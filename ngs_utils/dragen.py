import csv
import glob
import json
import re
from collections import defaultdict
import yaml
from os import listdir
from os.path import join, abspath, pardir, splitext, basename, dirname, realpath, isdir, isfile, exists

from ngs_utils.file_utils import adjust_path, verify_dir, file_exists, safe_mkdir, verify_file, add_suffix
from ngs_utils.logger import critical, debug, info, err, warn
from ngs_utils.Sample import BaseSample, BaseBatch, BaseProject
from natsort import natsort_keygen


class DragenSample(BaseSample):
    natsort_key = natsort_keygen()

    def __init__(self, **kwargs):
        BaseSample.__init__(self, **kwargs)  # name, dirpath, work_dir, bam, vcf, phenotype, normal_match
        self.qc_files = []
        self.bam = None


class DragenBatch(BaseBatch):
    def __init__(self, **kwargs):
        BaseBatch.__init__(self, **kwargs)
        self.somatic_caller = 'dragen'
        self.germline_caller = 'dragen'
        self.sv_caller = 'dragen'
        self.batch_qc_files = []
        self.somatic_vcf = join(self.parent_project.dir, self.name + '.hard-filtered.vcf.gz')
        self.sv_vcf = join(self.parent_project.dir, self.name + '.sv.vcf.gz')
        self.replay_file = join(self.parent_project.dir, self.name + '-replay.json')

    def find_somatic_vcf(self, silent=False):
        if isfile(self.somatic_vcf):
            verify_file(self.somatic_vcf, is_critical=True)
            if not silent:
                info(f'Found somatic VCF in <dragen-dir>/<tumor-name>-hard-filtered.vcf.gz: ' + self.somatic_vcf)

    def find_sv_vcf(self, silent=False):
        if isfile(self.sv_vcf):
            verify_file(self.sv_vcf, is_critical=True)
            if not silent:
                info(f'Found SV VCF in <dragen-dir>/<tumor-name>.sv.vcf.gz: ' + self.sv_vcf)

    def all_qc_files(self):
        return self.batch_qc_files + self.tumors[0].qc_files + self.normals[0].qc_files

    def add_tumor(self, name, rgid=None):
        sample = DragenSample(name=name, phenotype='tumor', batch=self, rgid=rgid)
        sample.bam = join(self.parent_project.dir, self.name + '_tumor.bam')
        self.tumors = [sample]
        if sample.name not in [s.name for s in self.parent_project.samples]:
            self.parent_project.samples.append(sample)
        return sample

    def add_normal(self, name, rgid=None):
        sample = DragenSample(name=name, phenotype='normal', batch=self, rgid=rgid)
        sample.bam = join(self.parent_project.dir, self.name + '.bam')
        self.normals = [sample]
        if sample.name not in [s.name for s in self.parent_project.samples]:
            self.parent_project.samples.append(sample)
        return sample


class DragenProject(BaseProject):
    @staticmethod
    def find_batches(input_dir, silent=False, include_samples=None, exclude_samples=None, parent_project=None):
        batch_by_name = dict()

        for replay_file in glob.glob(join(input_dir, '*-replay.json')):
            batch_name = basename(replay_file.split('-replay.json')[0])
            debug(f'Found somatic variants for batch {batch_name}')

            # Excluding/including based on batch name
            if exclude_samples and batch_name in exclude_samples:
                continue
            if include_samples and batch_name not in include_samples:
                continue

            # Reading *-replay.json to get the VCF sample names.
            # When DRAGEN is run with fastq inputs specificed directly with -1 -2 --tumor-fastq1 --tumor-fastq2,
            # user also must set the parameters --RGSM --RGSM-tumor, e.g.:
            # --RGSM P025_N --RGSM-tumor P025_T --output-directory /output/P025 --output-file-prefix P025
            # Those RGSM end up inside the output VCFs:
            # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  P025_N  P025_T
            # And the VCF files themselves, like all other output files, are prefixed with the --output-file-prefix
            # parameter e.g. P025.hard-filtered.vcf.gz and P025.sv.vcf.gz. The way to get the RGSM values
            # seems to be from the *-replay.json file:
            # {
            #   "dragen_config": [
            #     {
            #        "name": "RGSM",
            #        "value": "P025_N"
            #     },
            #     {
            #        "name": "RGSM-tumor",
            #        "value": "P025_T"
            #     },
            #   ]
            # }
            with open(replay_file) as f:
                replay_data = json.load(f)
            dragen_conf = {rec.get('name'): rec.get('value') for rec in replay_data['dragen_config']
                           if 'name' in rec and 'value' in rec}
            if 'RGSM' in dragen_conf:
                rgsm_n = dragen_conf['RGSM']
                rgsm_t = dragen_conf['RGSM-tumor']
                info(f'Normal RGSM as read from the replay file: {rgsm_n}')
                info(f'Tumour RGSM as read from the replay file: {rgsm_t}')
            else:
                # However if the FastQ input was in a form of a CSV file-lists files, e.g.:
                # --fastq-list inputs_normal.csv --tumor-fastq-list inputs_tumor.csv
                # then --RGSM(-tumor) will not be specificed and thus won't be reflected in "dragen_config", but instead
                # "dragen_config" will have:
                # {
                #   "dragen_config": [
                #     {
                #        "name": "tumor-fastq-list",
                #        "value": "\/data\/LIST\/fastqs_tumor.csv"
                #     },
                #     {
                #        "name": "fastq-list",
                #        "value": "\/data\/LIST\/fastqs_normal.csv"
                #     },
                #   ]
                # }
                # And the RGSM values will be present in the CSV files like the following:
                # RGID,RGSM,RGLB,Lane,Read1File,Read2File
                # TCTCTACT.CGCGGTTC.2,MDX190230_L1901041,UnknownLibrary,2,/data/FASTQ/191220_A00130_0127_BHMCNVDSXX/MDX190230_L1901041_S13_L002_R1_001.fastq.gz,/data/FASTQ/191220_A00130_0127_BHMCNVDSXX/MDX190230_L1901041_S13_L002_R2_001.fastq.gz
                #
                # Unfortunately those "\/data\/LIST\/fastqs_tumor.csv" paths are internal for a DRAGEN container and
                # don't correpond to a real file paths and even file names.
                #
                # Fortunately for the UMCCR ISL workflow, we know that those files are named as <output-prefix>-tumor.csv
                # and <output-prefix>-normal.csv:
                # https://github.com/umccr-illumina/stratus/blob/master/showcase/wfv.json#L276
                # And hopefully put into the DRAGEN output folder. Otherwise, we can't do much and will error out.
                #
                if 'tumor-fastq-list' not in dragen_conf or 'fastq-list' not in dragen_conf:
                    critical(f'Cannot find RGSM values: no RGSM(-tumor) nor (tumor-)fastq-list values '
                             f'found in {replay_file}')

                n_fastqs_fpath = join(input_dir, batch_name + '_normal.csv')
                t_fastqs_fpath = join(input_dir, batch_name + '_tumor.csv')
                if not verify_file(n_fastqs_fpath) or not verify_file(t_fastqs_fpath):
                    critical(f'Files {n_fastqs_fpath} or {t_fastqs_fpath} corresponding to (tumor-)fastq-list entries '
                             f'in {replay_file} are not found in the DRAGEN output folder {input_dir}. We expect the '
                             f'UMCCR ISL workflow to copy them there. If it did not happen, it means an issue with '
                             f'the workflow or that the run was not run with the UMCCR workflow at all. If the latter, '
                             f'please copy the fastq lists files manually into {input_dir} as <output-prefix>-tumor.csv>'
                             f'and <output-prefix>-normal.csv, or run with direct fastq inputs and RGSM specified:'
                             f' -1 -2 --tumor-fastq1 --tumor-fastq2 --RGSM --RGSM-tumor.')

                with open(n_fastqs_fpath) as f:
                    try:
                        rec = next(csv.DictReader(f))
                    except StopIteration:
                        critical(f'No lines in the fastqs input file {n_fastqs_fpath}; cannot read RGSM.')
                    else:
                        if 'RGSM' not in rec:
                            critical(f'Cannot find RGSM value in {n_fastqs_fpath}')
                        rgsm_n = rec['RGSM']
                with open(t_fastqs_fpath) as f:
                    try:
                        rec = next(csv.DictReader(f))
                    except StopIteration:
                        critical(f'No lines in the fastqs input file {n_fastqs_fpath}; cannot read RGSM.')
                    else:
                        if 'RGSM' not in rec:
                            critical(f'Cannot find RGSM value in {n_fastqs_fpath}')
                        rgsm_t = rec['RGSM']

                info(f'Normal RGSM: {rgsm_n} as read from {n_fastqs_fpath}')
                info(f'Tumour RGSM: {rgsm_t} as read from {t_fastqs_fpath}')

            # Excluding/including based on RG id
            if exclude_samples and (rgsm_t in exclude_samples or rgsm_n in exclude_samples):
                continue
            if include_samples and not (rgsm_t in include_samples or rgsm_n in include_samples):
                continue

            batch = DragenBatch(name=batch_name, parent_project=parent_project)
            batch_by_name[batch_name] = batch
            batch.add_tumor(batch_name, rgid=rgsm_t)
            batch.add_normal(batch_name + '_normal', rgid=rgsm_n)
            if exclude_samples and batch.normals[0].name in exclude_samples:
                continue
            batch.tumors[0].bam = join(input_dir, batch_name + '_tumor.bam')
            batch.normals[0].bam = join(input_dir, batch_name + '.bam')
            batch_by_name[batch_name] = batch

            batch.find_somatic_vcf(silent=silent)
            batch.find_germline_vcf(silent=silent)
            batch.find_sv_vcf(silent=silent)

            # populating qc files for multiqc:
            for suffix in [
                '.fragment_length_hist.csv',
                '.mapping_metrics.csv',
                '.ploidy_estimation_metrics.csv',
                '.time_metrics.csv',
                '.vc_metrics.csv',
                '.sv_metrics.csv',
            ]:
                qc_fpath = join(input_dir, f'{batch_name}{suffix}')
                if isfile(qc_fpath):
                    debug(f'Found QC file for batch {qc_fpath}')
                    batch.batch_qc_files.append(qc_fpath)
                else:
                    debug(f'Can\'t find QC file for batch {qc_fpath}')

            for suffix in [
                '.wgs_contig_mean_cov_{phenotype}.csv',
                '.wgs_coverage_metrics_{phenotype}.csv',
                '.wgs_fine_hist_{phenotype}.csv',
            ]:
                qc_fpath = join(input_dir, f'{batch_name}{suffix}')
                if isfile(qc_fpath.format(phenotype="normal")):
                    debug(f'Found QC file for normal sample {qc_fpath}')
                    batch.normals[0].qc_files.append(qc_fpath.format(phenotype="normal"))
                else:
                    debug(f'Can\'t find QC file for normal sample {qc_fpath}')
                if isfile(qc_fpath.format(phenotype="tumor")):
                    debug(f'Found QC file for tumor sample {qc_fpath}')
                    batch.tumors[0].qc_files.append(qc_fpath.format(phenotype="tumor"))
                else:
                    debug(f'Can\'t find QC file for tumor sample {qc_fpath}')

            debug(f'Found {len(batch.batch_qc_files)} batch QC files, '
                  f'{len(batch.tumors[0].qc_files)} tumor QC files, '
                  f'{len(batch.normals[0].qc_files)} normal QC files')
        return batch_by_name

    def __init__(self, input_dir=None, silent=False, include_samples=None, exclude_samples=None,
                 genome_build=None, **kwargs):
        BaseProject.__init__(self, input_dir=input_dir, **kwargs)
        self.genome_build = genome_build

        debug(f'Parsing project {input_dir}')
        self.batch_by_name = DragenProject.find_batches(self.dir, silent=silent,
            include_samples=include_samples, exclude_samples=exclude_samples, parent_project=self)

        if len(self.batch_by_name) == 1:
            self.project_name = list(self.batch_by_name.values())[0].name
        else:
            self.project_name = basename(input_dir)

    def add_batch(self, batch_name):
        batch = DragenBatch(name=batch_name, parent_project=self)
        self.batch_by_name[batch_name] = batch
        return batch

