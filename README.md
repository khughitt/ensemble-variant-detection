Ensemble Variant Detection
==========================

Andre Hennig, Keith Hughitt, Alexander Peltzer, Shrutii Sarda

Overview
--------

This project was started as part of a summer course in bioinformatics
and computational biology hosted by the [University of
Tübingen](http://www.ra.cs.uni-tuebingen.de/links/bioinformatik/welcome_e.html)
in collaboration with the [University of Maryland, College
Park](http://www.cbcb.umd.edu/) from August 4-9, 2014.

Requirements
------------

EVE is written in Python and requires a recent version of Python and several
Python libraries, as well as a number of command-line bioinformatics tools.

## Python

- [Python 3](https://www.python.org/downloads/)
- [NumPy](http://www.numpy.org/)
- [Pandas](http://pandas.pydata.org/)
- [PyVCF](http://pyvcf.readthedocs.org/en/latest/INTRO.html)
- [scikit-learn](http://scikit-learn.org/stable/)

## Bioinformatics tools

- [Burrows-Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/)
- [Genome Analysis Toolkit (GATK)](https://www.broadinstitute.org/gatk/)
- [SAMtools](http://samtools.sourceforge.net/)
- [Picard](http://picard.sourceforge.net/)

Installation
------------
@TODO

Configuration
-------------
@TODO

Usage
-----

## FASTQ input example

```
python eve.py -g path/to/annotations.gff    \
              -f path/to/genome.fasta       \
              -o output.vcf                 \
              reads_1.fastq reads_2.fastq
```

## BAM input example

```
python eve.py -g path/to/annotations.gff    \
              -f path/to/genome.fasta       \
              -o output.vcf                 \
              accepted_hits.bam
```

## Training example

```
python eve.py -g path/to/annotations.gff    \
              -f path/to/genome.fasta       \
              -o output.vcf                 \
              --train=actual_snps.vcf       \
              --num-threads=32              \
              reads_1.fastq reads_2.fastq
```

## A more complex example:

```
python eve.py --gff=path/to/annotations.gff            \
              --fasta=path/to/genome.fasta             \
              --mapper=bowtie2                         \
              --variant-detectors=gatk,mpileup,varscan \
              --working-directory=/scratch/eve         \
              --output=out.vcf                         \
              reads_1.fastq.gz reads_2.fastq.gz
```

TODO
----
- Add support for single-end reads
- Enable setting of Picard location
- Incorporate coverage,quality scores,sequence complexity and GC richness
  into classification.
- Check for FASTA indices
- unit testing / CI
- sphinx documentation
- setup.py

References
----------

-   Michael D Linderman, Tracy Brandt, Lisa Edelmann, Omar Jabado, Yumi
    Kasai, Ruth Kornreich, Milind Mahajan, Hardik Shah, Andrew
    Kasarskis, Eric E Schadt, (2014) Analytical Validation of Whole
    Exome And Whole Genome Sequencing For Clinical Applications. <em>Bmc
    Medical Genomics</em> <strong>7</strong> 20-NA
    <a href="http://dx.doi.org/10.1186/1755-8794-7-20">10.1186/1755-8794-7-20</a>
-   A. Talwalkar, J. Liptrap, J. Newcomb, C. Hartl, J. Terhorst, K.
    Curtis, M. Bresler, Y. S. Song, M. I. Jordan, D. Patterson, (2014)
    Smash: A Benchmarking Toolkit For Human Genome Variant Calling.
    <em>Bioinformatics</em>
    <a href="http://dx.doi.org/10.1093/bioinformatics/btu345">10.1093/bioinformatics/btu345</a>
