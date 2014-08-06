Ensemble Variant Detection
==========================

Andre Hennig, Keith Hughitt, Alexander Peltzer, Shrutii Sarda

Overview
--------

This project was started as part of a summer course in bioinformatics
and computational biology hosted by the university Tübingen from August
4-9, 2014.

Requirements
------------
- Python 3
- [PyYAML](http://pyyaml.org/)
- [Burrows-Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net/)

Installation
------------
@TODO

Configuration
-------------
@TODO

Usage
-----

Example usage:

```
python eve.py -i path/to/input/reads*.fastq \
              -g path/to/genome.gff         \
              -o output.vcf
```

A more complex example:

```
python eve.py --input=path/to/input/reads*.fastq       \
              --gff=path/to/genome.gff                 \
              --mapper=bowtie2                         \
              --variant-detectors=gatk,mpileup,varscan \
              --output=out.vcf
```

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
