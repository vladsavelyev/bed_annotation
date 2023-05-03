import sys
import math
from collections import OrderedDict
import os
from os.path import isfile, join, abspath, basename, dirname, getctime, getmtime, splitext, realpath
from subprocess import check_output

from ngs_utils.call_process import run
from ngs_utils import call_process
from ngs_utils.file_utils import intermediate_fname, iterate_file, splitext_plus, verify_file, adjust_path, add_suffix, \
    safe_mkdir, file_transaction, which, file_exists, open_gzipsafe, can_reuse
from ngs_utils.logger import info, critical, warn, err, debug
from ngs_utils import reference_data as ref
from ngs_utils.utils import md5


def get_chrom_order(genome=None, fai_fpath=None):
    chr_lengths = ref.get_chrom_lengths(genome, fai_fpath)
    chr_order = {c: i for i, (c, l) in enumerate(chr_lengths)}
    return chr_order


class SortableByChrom:
    def __init__(self, chrom, chrom_ref_order):
        self.chrom = chrom
        if isinstance(chrom_ref_order, dict):
            self.chrom_ref_order = chrom_ref_order.get(chrom, -1)
        else:
            self.chrom_ref_order = chrom_ref_order

    def get_key(self):
        return self.chrom_ref_order


class Region(SortableByChrom):
    def __init__(self, chrom, start, end, other_fields, chrom_ref_order):
        SortableByChrom.__init__(self, chrom, chrom_ref_order)
        self.start = int(start)
        self.end = int(end) if end is not None else None
        self.other_fields = list(other_fields)

    def get_key(self):
        return self.chrom_ref_order, self.start, self.end, self.other_fields

    def __str__(self):
        ts = [self.chrom, self.start, self.end] + self.other_fields
        return '"' + '\t'.join(map(str, ts)) + '"'

    def __repr__(self):
        return self.__str__()


def bam_to_bed(bam_fpath, to_gzip=True):
    debug('Converting the BAM to BED to save some memory.')  # from here: http://davetang.org/muse/2015/08/05/creating-a-coverage-plot-using-bedtools-and-r/
    bam_bed_fpath = splitext_plus(bam_fpath)[0] + ('.bed.gz' if to_gzip else '.bed')
    if can_reuse(bam_bed_fpath, bam_fpath):
        return bam_bed_fpath
    bedtools = which('bedtools')
    gzip = which('gzip')
    cmdline = '{bedtools} bamtobed -i {bam_fpath}'.format(**locals())
    cmdline += ' | {gzip}'.format(**locals()) if to_gzip else ''
    call_process.run(cmdline, output_fpath=bam_bed_fpath)
    return bam_bed_fpath


def bam_to_bed_nocnf(bam_fpath, bedtools='bedtools', gzip='gzip'):
    info('Converting the BAM to BED to save some memory.')  # from here: http://davetang.org/muse/2015/08/05/creating-a-coverage-plot-using-bedtools-and-r/
    bam_bed_fpath = splitext_plus(bam_fpath)[0] + '.bed.gz'
    cmdline = '{bedtools} bamtobed -i {bam_fpath} | {gzip} > {bam_bed_fpath}'.format(**locals())
    info(cmdline)
    os.system(cmdline)
    bam_bed_fpath = verify_file(bam_bed_fpath)
    if bam_bed_fpath:
        info('Done, saved to ' + bam_bed_fpath)
    else:
        err('Error, result is non-existent or empty')
    return bam_bed_fpath


def count_bed_cols(bed_fpath):
    with open(bed_fpath) as f:
        for l in f:
            if l and l.strip() and not l.startswith('#'):
                return len(l.split('\t'))
    # return len(next(dropwhile(lambda x: x.strip().startswith('#'), open(bed_fpath))).split('\t'))
    err('Empty bed file: ' + bed_fpath)
    return None


def merge_overlaps(work_dir, bed_fpath, distance=None):
    """Merge bed file intervals to avoid overlapping regions.
    Overlapping regions (1:1-100, 1:90-100) cause issues with callers like FreeBayes
    that don't collapse BEDs prior to using them.
    """
    output_fpath = intermediate_fname(work_dir, bed_fpath, 'merged')
    if isfile(output_fpath) and verify_file(output_fpath, cmp_f=bed_fpath):
        return output_fpath

    with file_transaction(work_dir, output_fpath) as tx:
        import pybedtools
        kwargs = dict(d=distance) if distance else dict()
        pybedtools.BedTool(bed_fpath).merge(**kwargs).saveas(tx)
    return output_fpath


def calc_region_number(bed_fpath):
    with open(bed_fpath) as f:
        return sum(1 for l in f if l.strip() and not l.strip().startswith('#'))


def get_genes_from_bed(bed_fpath, chrom_index=0, gene_index=3):
    gene_set = set()
    gene_list = list()
    with open(bed_fpath) as f:
        for line in f:
            if not line or not line.strip() or line.startswith('#'):
                continue

            tokens = line.replace('\n', '').split('\t')
            if len(tokens) <= gene_index or len(tokens) <= chrom_index:
                continue

            chrom = tokens[chrom_index]
            for gn in tokens[gene_index].split(','):
                if gn not in gene_set:
                    gene_set.add(gn)
                    gene_list.append(gn)

    return gene_set, gene_list


def cut(fpath, col_num, output_fpath=None):
    output_fpath = output_fpath or add_suffix(fpath, 'cut')
    if can_reuse(output_fpath, fpath):
        return output_fpath
    cmdline = 'cut -f' + ','.join(map(str, range(1, col_num + 1))) + ' ' + fpath + ' > ' + output_fpath
    call_process.run(cmdline, output_fpath=output_fpath)
    return output_fpath


def bgzip_and_tabix(fpath, reuse=False, tabix_parameters='', **kwargs):
    gzipped_fpath = join(fpath + '.gz')
    tbi_fpath = gzipped_fpath + '.tbi'

    if reuse and \
           file_exists(gzipped_fpath) and (getctime(gzipped_fpath) >= getctime(fpath) if file_exists(fpath) else True) and \
           file_exists(tbi_fpath) and getctime(tbi_fpath) >= getctime(gzipped_fpath):
        info('Actual compressed file and index exist, reusing')
        return gzipped_fpath

    info('Compressing and tabixing file, writing ' + gzipped_fpath + '(.tbi)')
    bgzip = which('bgzip')
    tabix = which('tabix')
    if not bgzip:
        err('Cannot index file because bgzip is not found')
    if not tabix:
        err('Cannot index file because tabix is not found')
    if not bgzip and not tabix:
        return fpath

    if isfile(gzipped_fpath):
        os.remove(gzipped_fpath)
    if isfile(tbi_fpath):
        os.remove(tbi_fpath)

    info('BGzipping ' + fpath)
    cmdline = '{bgzip} {fpath}'.format(**locals())
    call_process.run(cmdline)

    info('Tabixing ' + gzipped_fpath)
    cmdline = '{tabix} {tabix_parameters} {gzipped_fpath}'.format(**locals())
    call_process.run(cmdline)

    return gzipped_fpath


def filter_bed_with_gene_set(bed_fpath, gene_keys_set, output_fpath):
    with file_transaction(None, output_fpath) as tx:
        with open(bed_fpath) as inp, open(tx, 'w') as out:
            for l in inp:
                if l.strip('\n'):
                    chrom, start, end, gene = l.strip('\n').split('\t')
                    if (gene, chrom) in gene_keys_set:
                        out.write(l)


def sort_bed(input_bed_fpath, output_bed_fpath=None, work_dir=None, fai_fpath=None, chr_order=None, genome=None):
    input_bed_fpath = verify_bed(input_bed_fpath, is_critical=True)
    output_bed_fpath = adjust_path(output_bed_fpath) if output_bed_fpath \
        else intermediate_fname(work_dir, input_bed_fpath, 'sorted')

    debug('Sorting regions in ' + str(input_bed_fpath))
    if can_reuse(output_bed_fpath, input_bed_fpath):
        debug(output_bed_fpath + ' exists, reusing')
        return output_bed_fpath

    regions = []

    if not chr_order:
        if fai_fpath:
            fai_fpath = verify_file(fai_fpath)
        elif genome:
            fai_fpath = verify_file(ref.get_fai(genome))
        else:
            critical('Either of chr_order, fai_fpath, or genome build name must be specified')
        chr_order = get_chrom_order(fai_fpath=fai_fpath)

    with open(input_bed_fpath) as f:
        with file_transaction(work_dir, output_bed_fpath) as tx:
            with open(tx, 'w') as out:
                for l in f:
                    if not l.strip():
                        continue
                    if l.strip().startswith('#'):
                        out.write(l)
                        continue

                    fs = l.strip().split('\t')
                    chrom = fs[0]
                    start = int(fs[1])
                    end = int(fs[2])
                    other_fields = fs[3:]
                    order = chr_order.get(chrom, -1)
                    regions.append(Region(chrom, start, end, other_fields, order))

                for region in sorted(regions, key=lambda r: r.get_key()):
                    fs = [region.chrom, str(region.start), str(region.end)]
                    fs.extend(region.other_fields)
                    out.write('\t'.join(fs) + '\n')

    debug('Sorted ' + str(len(regions)) + ' regions, saved to ' + output_bed_fpath)
    return output_bed_fpath


def sort_bed_gsort(input_bed_fpath, output_bed_fpath=None, work_dir=None, fai_fpath=None, genome=None):
    input_bed_fpath = verify_bed(input_bed_fpath, is_critical=True)
    output_bed_fpath = adjust_path(output_bed_fpath) if output_bed_fpath \
        else intermediate_fname(work_dir, input_bed_fpath, 'sorted')

    debug('Sorting regions in ' + str(input_bed_fpath))
    if can_reuse(output_bed_fpath, input_bed_fpath):
        debug(output_bed_fpath + ' exists, reusing')
        return output_bed_fpath

    if fai_fpath:
        fai_fpath = verify_file(fai_fpath)
    elif genome:
        fai_fpath = verify_file(ref.get_fai(genome))
    else:
        critical('Either of fai_fpath or genome build name must be specified')

    with file_transaction(work_dir, output_bed_fpath) as tx:
        run('gsort {input_bed_fpath} {fai_fpath}'.format(**locals()), output_fpath=tx)

    return output_bed_fpath


def annotate_target(work_dir, target_bed, genome_build):
    output_fpath = intermediate_fname(work_dir, target_bed, 'ann')
    if can_reuse(output_fpath, target_bed):
        info(output_fpath + ' exists, reusing')
        return output_fpath

    bed_annotation = which('annotate_bed.py')
    if not bed_annotation:
        bed_annotation = which('bed_annotation')
        critical('Error: bed_annotation not found in PATH, please install `conda install -c vladsaveliev bed_annotation`.')

    cmdline = '{bed_annotation} {target_bed} -g {genome_build} -o {output_fpath}'.format(**locals())
    run(cmdline, output_fpath, stdout_to_outputfile=False)
    output_fpath = clean_bed(output_fpath, work_dir)
    return output_fpath


def calc_sum_of_regions(bed_fpath):
    total_bed_size = 0

    with open(bed_fpath) as f:
        for l in f:
            l = l.strip()
            if not l.startswith('#'):
                start, end = [int(f) for f in l.split('\t')[1:3]]
                total_bed_size += end - start

    return total_bed_size


def get_total_bed_size(bed_fpath, work_dir=None):
    import pybedtools
    if work_dir:
        pybedtools.set_tempdir(safe_mkdir(join(work_dir, 'pybedtools_tmp')))
    return sum(len(x) for x in pybedtools.BedTool(bed_fpath).merge())


def bedtools_version(bedtools):
    v = check_output([bedtools, '--version'])  # bedtools v2.24.0
    try:
        v = int(v.split(' ')[1].split('.')[1])
    except Exception:
        return None
    else:
        return v


this_py_fpath = splitext(__file__)[0] + '.py'
this_py_real_fpath = realpath(this_py_fpath)
project_dir = dirname(this_py_real_fpath)
code_base_path = abspath(join(abspath(project_dir)))


def get_padded_bed_file(work_dir, bed, padding, fai_fpath):
    genome_fpath = fai_fpath
    info('Making bed file for padded regions...')
    bedtools = which('bedtools')
    cmdline = '{bedtools} slop -i {bed} -g {genome_fpath} -b {padding}'.format(**locals())
    output_fpath = intermediate_fname(work_dir, bed, 'padded')
    call_process.run(cmdline, output_fpath=output_fpath)
    return output_fpath


def intersect_bed(work_dir, bed1, bed2):
    bed1_fname, _ = splitext_plus(basename(bed1))
    bed2_fname, _ = splitext_plus(basename(bed2))
    output_fpath = join(work_dir, bed1_fname + '__' + bed2_fname + '.bed')
    if can_reuse(output_fpath, [bed1, bed2]):
        return output_fpath
    bedtools = which('bedtools')
    cmdline = '{bedtools} intersect -u -a {bed1} -b {bed2}'.format(**locals())
    call_process.run(cmdline, output_fpath=output_fpath, checks=[call_process.file_exists_check])
    return output_fpath


def verify_bed(bed, description='', is_critical=False, silent=False):
    import pybedtools
    if isinstance(bed, pybedtools.BedTool):
        return bed

    fpath = adjust_path(bed)
    if not verify_file(fpath, description, is_critical=is_critical, silent=silent):
        return None

    error = BedFile(fpath).checkformat()
    if error:
        fn = critical if is_critical else err
        fn('Error: incorrect bed file format (' + fpath + '): ' + str(error) + '\n')
        return None

    return fpath


def clean_bed(bed_fpath, work_dir):
    clean_fpath = intermediate_fname(work_dir, bed_fpath, 'clean')

    if not can_reuse(clean_fpath, bed_fpath):
        import pybedtools
        pybedtools.set_tempdir(safe_mkdir(join(work_dir, 'pybedtools_tmp')))
        bed = pybedtools.BedTool(bed_fpath)
        bed = bed.filter(lambda x: x.chrom and
                         not any(x.chrom.startswith(e) for e in ['#', ' ', 'track', 'browser']))
        bed = bed.remove_invalid()
        with file_transaction(work_dir, clean_fpath) as tx_out_file:
            bed.saveas(tx_out_file)
        verify_bed(clean_fpath, is_critical=True)
        debug('Saved clean BED file into ' + clean_fpath)
    return clean_fpath


def check_md5(work_dir, fpath, file_type, silent=False):
    md5_fpath = join(work_dir, file_type + '_md5.txt')
    new_md5 = md5(fpath)
    info('md5 of ' + fpath + ' is ' + str(new_md5))
    prev_md5 = None
    if isfile(md5_fpath):
        with open(md5_fpath) as f:
            prev_md5 = f.read()
    else:
        info('Previous md5 file ' + md5_fpath + ' does not exist')
    info('Checking previous md5 from ' + md5_fpath + ': ' + str(prev_md5))

    if prev_md5 == new_md5:
        if not silent:
            debug('Reusing previous ' + file_type.upper() + ' files.')
        return True
    else:
        if not silent:
            info('Pre-processing input ' + file_type.upper() + ' file')
        if prev_md5:
            if not silent:
                info('Prev ' + file_type.upper() + ' md5: ' + str(prev_md5))
                info('New ' + file_type.upper() + ' md5: ' + str(new_md5))

        with open(md5_fpath, 'w') as f:
            f.write(str(new_md5))
        return False


class BedFile:
    def __init__(self, _filename):
        self.filename = _filename
        self.chrs = None
        self.nregions = None

    def checkformat(self):
        """************************************************************************************************************************************************************
        Task: checks the format of the bed file. The only requirements checked are that each line presents at least 3 tab separated columns, the
            two on the right must present integer values indicating the start/end position respectively. Right value must be greater than the
            left value.
        Outputs:
            err: string containing the detected error. Empty string in case of a correct format.
        ************************************************************************************************************************************************************"""

        fd = open_gzipsafe(self.filename)

        line = fd.readline()
        while line.startswith('#'):
            line = fd.readline()
        
        fields = line.split('\t')
        lc = 1
        error = ''

        # Checks that the two columns on the right contain integer values
        try:
            # Parses each line and checks that there are at least 3 fields, the two on the right containing integer values and being the right one
            # greater than the left one
            while line != '' and len(fields) > 2 and int(fields[1]) <= int(fields[2]):
                lc += 1
                line = fd.readline()
                fields = line.split('\t')
        except ValueError:
            error += 'Incorrect start/end values at line ' + str(lc) + '\n'
            error += 'Start/End coordinates must be indicated with integer values. The right value must be greater than the left value.\n'
            error += 'Line found: ' + line
            fd.close()

            return error

        # If it get to this point means that either the file ended or there is a line with less than 3 fields
        if line != '':
            error += 'Incorrect line format at line ' + str(lc) + '\n'
            error += 'At least three columns are expected in each line\n'
            error += 'The right value must be greater than the left value.\n'
            error += 'Line found: ' + line
            fd.close()

        return error

    def count_lines(self, filename=None):
        if filename is None:
            filename = self.filename
        return len(open(filename).readlines())


def calc_rate_within_normal(bases_by_depth, avg_depth, total_size):
    bases_within_normal = sum(
        bases
        for depth, bases
        in bases_by_depth.items()
        if math.fabs(avg_depth - depth) < 0.2 * avg_depth)

    return 1.0 * bases_within_normal / total_size if total_size else None


def calc_bases_within_threshs(bases_by_depth, total_size, depth_thresholds):
    bases_within_threshs = OrderedDict((depth, 0) for depth in depth_thresholds)
    rates_within_threshs = OrderedDict((depth, None) for depth in depth_thresholds)

    for depth, bases in bases_by_depth.items():
        for t in depth_thresholds:
            if depth >= t:
                bases_within_threshs[t] += bases
    for t in depth_thresholds:
        bs = bases_within_threshs[t]
        if total_size > 0:
            rate = 1.0 * bases_within_threshs[t] / total_size
            if rate > 1:
                critical('Error: rate is > 1: rate = ' + str(rate) + ', bases = ' + str(bs) + ', size = ' + str(total_size))
            rates_within_threshs[t] = rate

    return bases_within_threshs, rates_within_threshs
