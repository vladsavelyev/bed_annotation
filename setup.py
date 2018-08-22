#!/usr/bin/env python
import sys
py_v = sys.version_info[:2]
if py_v < (3, 6):
    sys.exit('Only Python 3.6 and higher are supported. Current version: ' + '.'.join(py_v))

from os.path import join

name = 'bed_annotation'
script_name = 'bed_annotation'
package_name = 'bed_annotation'

version = os.environ.get('TRAVIS_TAG', 'dev')

with open('requirements.txt') as f:
    reqs = f.read().strip().split('\n')

from setuptools import setup, find_packages
setup(
    name=name,
    version=version,
    author='Vlad Saveliev and Alla Mikheenko',
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
            'ensembl/hg19/canon_transcripts_hg19_ensembl.txt',
            'ensembl/hg38/ensembl.bed.gz',
            'ensembl/hg38/ensembl.bed.gz.tbi',
            'ensembl/hg38/canon_transcripts_hg38_ensembl.txt',
            'ensembl/mm10/ensembl.bed.gz',
            'ensembl/mm10/ensembl.bed.gz.tbi',
            'ensembl/canon_cancer_replacement.txt',
        ]
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
    install_requires=reqs,
    setup_requires=[],
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
