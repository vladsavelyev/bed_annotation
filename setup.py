#!/usr/bin/env python
from os.path import join
import sys
from setuptools import setup, find_packages

py_v = sys.version_info[:2]
if py_v < (3, 6) or py_v > (3, 10):
    sys.exit('Python versions 3.6 through 3.10 are supported. Current version: ' + '.'.join(py_v))

name = script_name = package_name = 'bed_annotation'

setup(
    name=name,
    version='1.2.0',
    author='Vlad Savelyev and Alla Mikheenko',
    author_email='vladislav.sav@gmail.com',
    description='Genome capture target coverage evaluation tool',
    long_description=(open('README.md').read()),
    keywords='bioinformatics',
    url='https://github.com/vladsaveliev/bed_annotation',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        package_name: [
            'ensembl/hg19/ensembl.bed.gz',
            'ensembl/hg19/ensembl.bed.gz.tbi',
            'ensembl/hg19/appris_data.principal.txt',
            'ensembl/hg38/ensembl.bed.gz',
            'ensembl/hg38/ensembl.bed.gz.tbi',
            'ensembl/hg38/appris_data.principal.txt',
            'ensembl/mm10/ensembl.bed.gz',
            'ensembl/mm10/ensembl.bed.gz.tbi',
        ],
        'ngs_utils': [
            'reference_data/fai/GRCh37.fa.fai',
            'reference_data/fai/hg19.fa.fai',
            'reference_data/fai/hg19-noalt.fa.fai',
            'reference_data/fai/hg19-chr21.fa.fai',
            'reference_data/fai/hg38.fa.fai',
            'reference_data/fai/hg38-noalt.fa.fai',
            'reference_data/fai/mm10.fa.fai',
        ],
    },
    include_package_data=True,
    zip_safe=False,
    scripts=[
        join('scripts', script_name),
        join('scripts', 'extract_features'),
        join('scripts', 'generate_cds_bed.py'),
        join('scripts', 'generate_ensembl_data.py'),
        join('scripts', 'generate_refseq_data.py'),
    ],
    install_requires=[
        'pip',
        'setuptools>=18.5',
        'pybedtools',
        'cython',
        'numpy',
        'joblib',
        'gffutils',
        'click',
    ],
    classifiers=[
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
)
