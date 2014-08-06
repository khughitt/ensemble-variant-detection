#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import glob
import argparse
from eve import detectors

class EVE(object):
    """Ensemble Variant Detection"""
    def __init__(self, argv):
        self.args = self.parse_args(argv)

    def run(self):
        """Main application process"""
        # load mappers

        # load detectors
        self.detectors = []

        for detector in self.args.variant_detectors.split(','):
            conf = os.path.join('config', '%s.yaml' % detector)
            self.detectors.append(detectors.Detector(conf))

    def parse_args(self, argv):
        """Parses input arguments"""
        parser = argparse.ArgumentParser(
                description='Ensemble Variant Detection')
        parser.add_argument('-i', '--input', required=True,
                            help=('Wildcard string specifying location of '
                                  'input FASTQ files to use'))
        parser.add_argument('-g', '--gff', required=True,
                            help='Location of GFF annotation file to use.')
        parser.add_argument('-m', '--mapper', default='bwa',
                            help='Mapper to use for read alignment')
        parser.add_argument('-d', '--variant-detectors',
                            default='gatk,mpileup,varscan',
                            help=('Comma-separated list of the variant '
                                  'detectors to be used.'))
        args = parser.parse_args()

        # validate input arguments
        args.input_list = glob.glob(args.input)
        if len(args.input_list) is 0:
            raise IOError("Invalid input filepath string specified")
        if not os.path.isfile(args.gff):
            raise IOError("Invalid GFF filepath specified")

        return args

if __name__ == '__main__':
    app = EVE(sys.argv)
    sys.exit(app.run())
