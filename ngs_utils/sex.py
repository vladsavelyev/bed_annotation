from __future__ import division
import pybedtools
from os.path import join, dirname, abspath
from pybedtools import BedTool

from ngs_utils.bed_utils import filter_bed_with_gene_set, sort_bed, get_total_bed_size
from ngs_utils.logger import err, info, debug, critical, warn
from ngs_utils.file_utils import verify_file, safe_mkdir, file_transaction, can_reuse
from ngs_utils.parallel import parallel_view, ParallelCfg
from ngs_utils.reporting.reporting import get_val, get_float_val, get_int_val
from ngs_utils.sambamba import index_bam, sambamba_depth
import ngs_utils.reference_data as ref


def get_gender(genome, bam_fpath, bed_fpath, sample, avg_depth):
    gender = None
    chrom_lengths = ref.get_chrom_lengths(genome)
    chrom_names = [chrom for chrom, length in chrom_lengths]
    if 'Y' in chrom_names or 'chrY' in chrom_names:
        gender = determine_sex(sample.work_dir, bam_fpath, avg_depth, genome, bed_fpath)
        if gender:
            with open(join(safe_mkdir(sample.dirpath), 'gender.txt'), 'w') as f:
                f.write(gender[0].upper())
    return gender


chry_key_regions_by_genome = {
    'hg19': join(dirname(abspath(__file__)), 'gender', 'chrY.hg19.bed'),
    'hg38': join(dirname(abspath(__file__)), 'gender', 'chrY.hg38.bed'),
}
AVG_DEPTH_THRESHOLD_TO_DETERMINE_SEX = 5
FEMALE_Y_COVERAGE_FACTOR = 10.0


def determine_sex(work_dir, bam_fpath, avg_depth, genome, target_bed=None):
    debug()
    debug('Determining sex')
    pybedtools.set_tempdir(safe_mkdir(join(work_dir, 'pybedtools_tmp')))

    male_bed = None
    for k in chry_key_regions_by_genome:
        if k in genome:
            male_bed = BedTool(chry_key_regions_by_genome.get(k))
            break
    if not male_bed:
        warn('Warning: no male key regions for ' + genome + ', cannot identify sex')
        return None

    male_area_size = get_total_bed_size(male_bed)
    debug('Male region total size: ' + str(male_area_size))

    if target_bed:
        target_male_bed = join(work_dir, 'male.bed')
        with file_transaction(work_dir, target_male_bed) as tx:
            BedTool(target_bed).intersect(male_bed).merge().saveas(tx)
        target_male_area_size = get_total_bed_size(target_male_bed)
        if target_male_area_size == 0:
            debug('The male non-PAR region does not overlap with the capture target - cannot determine sex.')
            return None
        male_bed = target_male_bed
    else:
        debug('WGS, determining sex based on chrY key regions coverage.')

    info('Detecting sex by comparing the Y chromosome key regions coverage and average coverage depth.')
    if not bam_fpath:
        critical('BAM file is required.')
    index_bam(bam_fpath)

    chry_mean_coverage = _calc_mean_coverage(work_dir, male_bed, bam_fpath, 1)
    debug('Y key regions average depth: ' + str(chry_mean_coverage))
    avg_depth = float(avg_depth)
    debug('Sample average depth: ' + str(avg_depth))
    if avg_depth < AVG_DEPTH_THRESHOLD_TO_DETERMINE_SEX:
        debug('Sample average depth is too low (less than ' + str(AVG_DEPTH_THRESHOLD_TO_DETERMINE_SEX) +
             ') - cannot determine sex')
        return None

    if chry_mean_coverage == 0:
        debug('Y depth is 0 - it\s female')
        sex = 'F'
    else:
        factor = avg_depth / chry_mean_coverage
        debug('Sample depth / Y depth = ' + str(factor))
        if factor > FEMALE_Y_COVERAGE_FACTOR:  # if mean target coverage much higher than chrY coverage
            debug('Sample depth is more than ' + str(FEMALE_Y_COVERAGE_FACTOR) + ' times higher than Y depth - it\s female')
            sex = 'F'
        else:
            debug('Sample depth is not more than ' + str(FEMALE_Y_COVERAGE_FACTOR) + ' times higher than Y depth - it\s male')
            sex = 'M'
    debug('Sex is ' + sex)
    debug()
    return sex


def _calc_mean_coverage(work_dir, bed_file, bam_fpath, threads):
    sambamba_depth_file = sambamba_depth(work_dir, bed_file, bam_fpath, threads=threads)

    mean_cov = []
    mean_cov_col = None
    total_len = 0
    with open(sambamba_depth_file) as bedcov_file:
        for line in bedcov_file:
            if line.startswith('#'):
                mean_cov_col = line.split('\t').index('meanCoverage')
                continue
            line_tokens = line.replace('\n', '').split()
            start, end = map(int, line_tokens[1:3])
            size = end - start
            mean_cov.append(float(line_tokens[mean_cov_col]) * size)
            total_len += size
    mean_cov = sum(mean_cov) / total_len if total_len > 0 else 0
    return mean_cov
