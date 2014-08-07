#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ensemble Variant Detection (EVE)

This file contains the main class and execution logic for the EVE variant
detection pipeline.
"""
import os
import sys
import glob
import logging
import argparse
import datetime
import platform
import subprocess
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
        self.log_system_info()

        # check for FASTA index and create if necessary
        self.check_fasta_index()

        # load mapper
        if 'bam' not in self.args:
            # split fastq reads into two variables
            (reads1, reads2) = self.args.input_reads

            # output file
            prefix = os.path.basename(
                        os.path.commonprefix([reads1, reads2])).strip("_")
            filename = "aln_%s.sam" % prefix
            outfile = os.path.join(self.working_dir, 'mapped', filename)

            self.mapper = mappers.BWAMemMapper(self.args.fasta, reads1, reads2,
                                               outfile, self.args.max_threads)

    def run(self):
        """Main application process"""
        # map reads
        if 'bam' not in self.args:
            logging.info("Mapping reads")
            self.args.bam = self.mapper.run()

        # load detectors
        self.load_detectors()

        # run variant detectors
        logging.info("Running variant detection algorithms")

        # TESTING (GATK)
        logging.info("Running GATK")
        self.detectors[0].run()

        # normalize output from variant detectors and read in as either a NumPy
        # matrix or pandas DataFrame

        # run classifier
        # (http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)

        # output final VCF

    def check_fasta_index(self):
        """Checks for a valid FASTA index and creates one if needed"""
        if not os.path.exists("%s.fai" % self.args.fasta):
            # FASTA indexing command
            cmd = "samtools faidx %s" % self.args.fasta

            logging.info("Creating a FASTA index")
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)

    def load_detectors(self):
        """Loads the variant detector instances"""
        # available detectors
        mapping = {
            'mpileup': detectors.MpileupDetector,
            'gatk'   : detectors.GATKDetector,
            'varscan': detectors.VarScanDetector
        }

        # load detectors
        self.detectors = []

        for detector in self.args.variant_detectors.split(','):
            conf = os.path.join('config', 'detectors', '%s.cmd' % detector)
            cls = mapping[detector]

            self.detectors.append(cls(
                self.args.bam, self.args.fasta, conf, self.working_dir,
                self.args.max_threads
            ))

    def create_working_directories(self):
        """Creates directories to output intermediate files into"""
        now = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')

        self.working_dir = os.path.join(self.args.working_directory, now)

        for subdir in ['mapped', 'vcf', 'mpileup']:
            path = os.path.join(self.working_dir, subdir)
            if not os.path.isdir(path):
                os.makedirs(path)

    def initialize_logger(self):
        """Initializes a logger instance"""
        logging.basicConfig(level=logging.DEBUG,
                format='(%(asctime)s)[%(levelname)s] %(message)s',
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

    def log_system_info(self):
        """Prints system information to the log.
           Code adapted from the SunPy project."""
        # EVE version
        from eve import __version__ as eve_version
        logging.info("Starting EVE %s" % eve_version)

        # Time
        now = datetime.datetime.utcnow().strftime("%A, %d. %B %Y %I:%M%p UT")
        logging.info("Time: %s" % now)

        # Platform
        system = platform.system()
        processor = platform.processor()

        if system == "Linux":
            distro = " ".join(platform.linux_distribution())
            logging.debug("OS: %s (Linux %s %s)" %  (
                distro, platform.release(), processor))
        elif system == "Darwin":
            logging.debug("OS: Mac OS X %s (%s)" % (
                platform.mac_ver()[0], processor)
            )
        elif system == "Windows":
            logging.debug("OS: Windows %s %s (%s)" %  (
                platform.release(), platform.version(), processor))
        else:
            logging.debug ("Unknown OS (%s)" % processor)

        # Architecture
        logging.debug('Architecture: %s' % platform.architecture()[0])

        # Python version
        logging.debug("Python %s" % platform.python_version())

        # Check python dependencies
        try:
            from numpy import __version__ as numpy_version
        except ImportError:
            numpy_version = "NOT INSTALLED"

        try:
            from scipy import __version__ as scipy_version
        except ImportError:
            scipy_version = "NOT INSTALLED"

        try:
            from sklearn import __version__ as sklearn_version
        except ImportError:
            sklearn_version = "NOT INSTALLED"

        logging.debug("NumPy: %s" % numpy_version)
        logging.debug("SciPy: %s" % scipy_version)
        logging.debug("Scikit-Learn: %s" % sklearn_version)

        # @TODO: command-line tool versions (SAMtools, etc)

    def parse_args(self, argv):
        """Parses input arguments"""
        parser = argparse.ArgumentParser(
                description='Ensemble Variant Detection')
        parser.add_argument('input_reads', nargs='+',
                            help=('Input paired-end Illumina reads or '
                                  'alignment. Supported file formats include '
                                  '.fastq, .fastq.gz, and .bam'))
        parser.add_argument('-o', '--output',
                            help='Location to save final VCF output to.')
        parser.add_argument('-f', '--fasta', required=True,
                            help='Location of genome sequence file to use.')
        parser.add_argument('-g', '--gff', required=True,
                            help='Location of GFF annotation file to use.')
        parser.add_argument('-m', '--mapper', default='bwa',
                            help='Mapper to use for read alignment')
        parser.add_argument('-t', '--max-threads', default='4',
                            help='Maximum number of threads to use')
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
            args.bam = args.input_reads[0]

        return args

if __name__ == '__main__':
    app = EVE(sys.argv)
    sys.exit(app.run())
