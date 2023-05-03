def is_autosomal(chrom):
    """Keep chromosomes that are a digit 1-22, or chr prefixed digit chr1-chr22
    """
    try:
        int(chrom)
        return True
    except ValueError:
        try:
            int(str(chrom.replace("chr", "")))
            return True
        except ValueError:
            return False

def is_sex(chrom):
    return chrom in ["X", "chrX", "Y", "chrY"]

def is_mitochondrial(chrom):
    return chrom in ["MT", "chrM", "chrMT"]

def is_autosomal_or_x(chrom):
    return is_autosomal(chrom) or chrom in ["X", "chrX"]

def is_autosomal_or_sex(chrom):
    return is_autosomal_or_x(chrom) or is_sex(chrom)

def is_alt(chrom):
    """Check that a chromosome is on 1-22, X, Y, MT.
    """
    return not is_autosomal_or_sex(chrom) and not is_mitochondrial(chrom)
