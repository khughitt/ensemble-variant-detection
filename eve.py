#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import logging
import argparse
import datetime
from eve import detectors,mappers

class EVE(object):
    """Ensemble Variant Detection"""
    def __init__(self, argv):
        # parse arguments
        self.args = self.parse_args(argv)

        # create working directories
        self.create_working_directories()

        # initialize logger
        self.initialize_logger()

        # load mapper
        if args.input_type == 'fastq':
            self.mapper = mappers.BWAMapper()

        # load detectors
        self.detectors = []

        for detector in self.args.variant_detectors.split(','):
            conf = os.path.join('config', 'detectors', '%s.yaml' % detector)
            self.detectors.append(detectors.Detector(conf))

    def create_working_directories(self):
        """Creates directories to output intermediate files into"""
        now = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')

        self.working_dir = os.path.join(self.args.working_directory, now)

        for subdir in ['mapped', 'vcf']:
            path = os.path.join(self.working_dir, subdir)
            if not os.path.isdir(path):
                os.makedirs(path)

    def initialize_logger(self):
        """Initializes a logger instance"""
        logging.basicConfig(level=logging.INFO,
                format='%(asctime)s [%(levelname)s] %(message)s',
                filename=os.path.join(self.working_dir, 'eve.log'))

        # log to console as well
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)

        # set a format which is simpler for console use
        formatter = logging.Formatter('[%(levelname)s] %(message)s')

        # tell the handler to use this format
        console.setFormatter(formatter)

        # add the handler to the root logger
        logging.getLogger('').addHandler(console)

    def run(self):
        """Main application process"""
        # map reads
        if args.input_type == 'fastq':
            logging.info("Mapping reads")
            self.mapper.run(self.args.input_reads)

        # run variant detectors

        # normalize output from variant detectors and read in as either a NumPy
        # matrix or pandas DataFrame

        # run classifier
        # (http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)

        # output final VCF


    def parse_args(self, argv):
        """Parses input arguments"""
        parser = argparse.ArgumentParser(
                description='Ensemble Variant Detection')
        parser.add_argument('input_reads', nargs='+',
                            help=('Input paired-end Illumina reads or '
                                  'alignment. Supported file formats include '
                                  '.fastq, .fastq.gz, and .bam')
        parser.add_argument('-o', '--output',
                            help='Location to save final VCF output to.')
        parser.add_argument('-g', '--gff', required=True,
                            help='Location of GFF annotation file to use.')
        parser.add_argument('-m', '--mapper', default='bwa',
                            help='Mapper to use for read alignment')
        parser.add_argument('-d', '--variant-detectors',
                            default='gatk,mpileup,varscan',
                            help=('Comma-separated list of the variant '
                                  'detectors to be used.'))
        parser.add_argument('-w', '--working-directory',
                            default='output/',
                            help='Location to store intermediate files')
        args = parser.parse_args()

        # validate input arguments
        if len(args.input_reads) > 2:
            raise IOError("Too many input arguments specified")
        if not os.path.isfile(args.gff):
            raise IOError("Invalid GFF filepath specified")

        # determine input type (FASTQ or BAM)
        if len(args.input_reads) == 1 and args.input_reads[0].endswith('.bam'):
            args.input_type = 'bam'
        else:
            args.input_type = 'fastq'

        return args

if __name__ == '__main__':
    app = EVE(sys.argv)
    sys.exit(app.run())
