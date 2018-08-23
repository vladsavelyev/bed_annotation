# BED Annotation

[![Build Status](https://travis-ci.org/vladsaveliev/bed_annotation.svg?branch=master)](https://travis-ci.org/vladsaveliev/bed_annotation)
[![Anaconda-Server Badge](https://anaconda.org/vladsaveliev/bed_annotation/badges/installer/conda.svg)](https://conda.anaconda.org/vladsaveliev)

A tool that assigns gene names to regions in a BED file based on Ensembl genomic features overlap.

### Installation

```
conda install -c vladsaveliev bed_annotation
```

### Usage

```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed
``` 

The script checks each BED region against the Ensembl genomic features database, and writes a BED file in a standardized format with a gene symbol, strand and exon rank in 4-6th columns:

`INPUT.bed`:

```
chr1    69090   70008
chr1    367658  368597
```

`OUTPUT.bed`:

```
chr1    69090   70008   OR4F5   1       +
chr1    367658  368597  OR4F29  1       +
```

Available genomes (to provide with `-g`): GRCh37, hg19, hg38.

#### Transcripts order

The piority for choosing transcripts for annotation is the following:
- Overlap % with transcript
- Overlap % with CDS
- Overlap % with exons
- Biotype (`protein_coding` > others > `*RNA` > `*_decay` > `sense_*` > `antisense` > `translated_*` > `transcribed_*`)
- TSL (1 > NA > others > 2 > 3 > 4 > 5)
- Presence of a HUGO gene symbol
- Is cancer canonical
- Transcript size

#### Extended annotation

Use `--extended` option to report extra columns with details on features, biotype, overlapping transcripts and overlap sizes:

```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed --extended
```

`OUTPUT.bed`:

```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon  Strand  Feature Biotype         Ensembl_ID      TSL HUGO    Tx_overlap_% Exon_overlaps_% CDS_overlaps_% Ori_Fields
chr1    69090   70008   OR4F5   1     +       capture protein_coding  ENST00000335137 NA  OR4F5   100.0        100.0           99.7
chr1    367658  368597  OR4F29  1     +       capture protein_coding  ENST00000426406 NA  OR4F29  100.0        100.0           99.7
```

#### Ambuguous annotations

Regions may overlap mltiple genes. The `--ambiguities` controls how the script resolves such ambiguities

- `--ambiguities all` -- report all reliable overlaps (in order in the "priority" section, see above)
- `--ambiguities all_ask` -- stop execution and ask user which annotation to pick
- `--ambiguities best_all` (default) -- find the best overlap, and if there are several equally good, report all (in terms of the "priority" above)
- `--ambiguities best_ask` -- find the best overlap, and if there are several equally good, ask user
- `--ambiguities best_one` -- find the best overlap, and if there are several equally good, report any of them

Note that the first 4 options might output multiple lines per region, e.g.:

```
annotate_bed.py INPUT.bed -g hg19 -o OUTPUT.bed --extended --ambiguities best_all
```

`OUTPUT.bed`:

```
## Tx_overlap_%: part of region overlapping with transcripts
## Exon_overlaps_%: part of region overlapping with exons
## CDS_overlaps_%: part of region overlapping with protein coding regions
#Chrom  Start   End     Gene    Exon    Strand  Feature Biotype Ensembl_ID      TSL     HUGO    Tx_overlap_%    Exon_overlaps_% CDS_overlaps_%
chr1    69090   70008   OR4F5   1       +       capture protein_coding  ENST00000335137 NA      OR4F5   100.0   100.0   100.0
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000426406 NA      OR4F29  100.0   100.0   100.0
chr1    367658  368597  OR4F29  1       +       capture protein_coding  ENST00000412321 NA      OR4F29  100.0   100.0   100.0
```

#### Other options

- `--coding-only`: take only the features of type `protein_coding` for annotation
- `--high-confidence`: annotate with only high confidence regions (TSL is 1 or NA, with HUGO symbol, total overlap size > 50%)
- `--canonical`: use only canonical transcripts to annotate (which to the most part means the longest transcript, by SnpEff definition)
- `--short`: add only the 4th "Gene" column (outputa 4-col BED file instead of 6-col)
- `--output-features`: good for debugging. Under each BED file region, also output Ensemble featues that were used to annotate it
