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
        self.args = self.parse_args(argv)

        # create working directories
        self.create_working_directories()

        # initialize logger
        self.initialize_logger()

        # load mapper
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
            if not os.path.isdir(subdir):
                os.makedirs(subdir)

    def initialize_logger(self):
        """Initializes a logger instance"""
        logging.basicConfig(level=logging.INFO,
                            filename=os.path.join(self.working_dir, 'eve.log')

        # log to console as well
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(formatter)

        # add the handler to the root logger
        logging.getLogger('').addHandler(console)

    def run(self):
        """Main application process"""
        # map reads
        logging.info("Mapping reads")

        self.mapper.run(args.input_reads)

    def parse_args(self, argv):
        """Parses input arguments"""
        parser = argparse.ArgumentParser(
                description='Ensemble Variant Detection')
        parser.add_argument('input-reads', nargs='+',
                            help='Input paired-end Illumina reads')
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
                            default='output',
                            help='Location to store intermediate files')
        args = parser.parse_args()

        # validate input arguments
        if len(args.input_reads) > 2:
            raise IOError("Too many input arguments specified")
        if not os.path.isfile(args.gff):
            raise IOError("Invalid GFF filepath specified")

        return args

if __name__ == '__main__':
    app = EVE(sys.argv)
    sys.exit(app.run())
