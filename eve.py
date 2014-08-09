#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ensemble Variant Detection (EVE)

This file contains the main class and execution logic for the EVE variant
detection pipeline.
"""
import os
import sys
import vcf
import glob
import logging
import argparse
import datetime
import platform
import subprocess
import configparser
from pandas import DataFrame
from eve import detectors,mappers


class EVE(object):
    """Ensemble Variant Detection"""
    def __init__(self, argv):
        # parse arguments
        self.args = self.parse_args(argv)

        # determine working/output directory to use
        if self.args.output_dir == 'output/{timestamp}'
            now = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
            self.output_dir = os.path.join(self.args.output_dir, now)
        else:
            self.output_dir = self.args.output_dir

        # create working directories
        self.create_output_directories()

        # initialize logger
        self.initialize_logger()
        self.log_system_info()

        # load configuration
        self.load_config()

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
            outfile = os.path.join(self.output_dir, 'mapped', filename)

            self.mapper = mappers.BWAMemMapper(self.args.fasta, reads1, reads2,
                                               outfile, self.args.num_threads)

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

        # @TODO: parallelize
        vcf_files = []

        for detector in self.detectors:
            output = detector.run()

            if isinstance(output, list):
                vcf_files += output
            else:
                vcf_files.append(output)

        # normalize output from variant detectors and read in as either a NumPy
        # matrix or pandas DataFrame
        df = self.combine_vcfs(vcf_files)
        df.to_csv(os.path.join(self.output_dir, "combined.csv"),
                  index_label='position')

        # run classifier
        # (http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html)
        if self.args.training_set:
            self.build_training_set(df)

        # output final VCF

    def build_training_set(self, df):
        """Adds actual values to the end of the combined dataset"""
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier

        # load VCF containing true answers
        # For now, assuming Genome in a Bottle VCF...
        reader = vcf.Reader(open(self.args.training_set))

        # Add "truth" values
        df['actual'] = np.repeat(float('nan'), df.shape[0])

        for record in reader:
            #if record.POS in df['position'].values:
            if record.POS in df.index:
                # there is probably a cleaner way to do this, but I can't
                # think of it right now...
                #df.loc[df[df.position == pos].index, 'actual'] = record.ALT[0]
                df.loc[df[df.index == record.POS].index, 'actual'] = record.ALT[0]

        df.to_csv(os.path.join(self.output_dir, "combined_training_set.csv"),
                  index_label='position')

    def combine_vcfs(self, vcf_files):
        """Parses a collection of VCF files and creates a single matrix
        containing the calls for each position observed by any of the detection
        algorithms."""
        logging.info("Combining output from variant detection tools")

        # Load VCF files
        unfiltered_vcf_readers = {}

        for filename in vcf_files:
            name = os.path.splitext(os.path.basename(filename))[0]

            # if indels, skip...
            if name == 'varscan_indels':
                continue

            unfiltered_vcf_readers[name] = vcf.Reader(open(filename))

        # Remove filtered entries and (for now) non-SNPs
        filtered_vcf_readers = {}

        for name, reader in unfiltered_vcf_readers.items():
            filtered = []

            for record in reader:
                if record.is_snp:
                    if not record.FILTER or len(record.FILTER) == 0:
                        filtered.append(record)

            filtered_vcf_readers[name] = filtered

        # Get calls for each algorithm at each position observed by any of the
        # tools
        combined_dict = {name:{} for name in filtered_vcf_readers.keys()}
        combined_dict['depth'] = {}

        for name,reader in filtered_vcf_readers.items():
            # create a new column for the quality scores
            qual_name = name + "_qual"
            combined_dict[qual_name] = {}

            for record in reader:
                # @TODO: decide how to deal with multiple alleles
                # i.e.: len(record.ALT) > 1

                # Determine quality score to use
                try:
                    # GATK
                    qual_score = record.INFO['QD']
                except KeyError:
                    try:
                        # VarScan
                        # http://varscan.sourceforge.net/support-faq.html#output-confidence
                        qual_score = record.samples[0]['GQ']
                    except:
                        # mpileup
                        # Also contains the Genotype Quality score used for
                        # VarScan above...
                        qual_score = record.QUAL / record.INFO['DP']

                # Determine read depth (only need once...)
                try:
                    # GATK / mpileup
                    depth = record.INFO['DP']
                except KeyError:
                    # VarScan
                    depth = record.INFO['ADP']

                # Add entries to the dictionary
                combined_dict['depth'][record.POS] = depth
                combined_dict[name][record.POS] = record.ALT[0]
                combined_dict[qual_name][record.POS] = qual_score

        # Convert to a DataFrame
        return DataFrame.from_dict(combined_dict)

    def check_fasta_index(self):
        """Checks for a valid FASTA index and creates one if needed"""
        # FAI index
        if not os.path.exists("%s.fai" % self.args.fasta):
            # FASTA indexing command
            cmd = "samtools faidx %s" % self.args.fasta

            logging.info("Creating a FASTA index")
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)

        # Sequence dictionary
        base_filename = os.path.splitext(self.args.fasta)[0]
        seqdict = "%s.dict" % base_filename

        if not os.path.exists(seqdict):
            cmd = "java -jar CreateSequenceDictionary.jar R=%s O=%s" % (
                self.args.fasta, seqdict
            )
            logging.info("Creating a FASTA sequence dictionary")
            logging.debug(cmd)
            subprocess.call(cmd, shell=True)

    def load_detectors(self):
        """Loads the variant detector instances"""
        # available detectors
        mapping = {
            'mpileup': (detectors.MpileupDetector, None),
            'gatk'   : (detectors.GATKDetector,
                        self.config['tools']['GATK_jar']),
            'varscan': (detectors.VarScanDetector,
                        self.config['tools']['VarScan_jar'])
        }

        # load detectors
        self.detectors = []

        for detector in self.args.variant_detectors.split(','):
            conf = os.path.join('config', 'detectors', '%s.cmd' % detector)
            (cls, location) = mapping[detector]

            self.detectors.append(cls(
                self.args.bam, self.args.fasta, conf, self.output_dir,
                self.args.num_threads, location
            ))

    def create_output_directories(self):
        """Creates directories to output intermediate files into"""

        for subdir in ['mapped', 'vcf', 'mpileup']:
            path = os.path.join(self.output_dir, subdir)
            if not os.path.isdir(path):
                os.makedirs(path)

    def initialize_logger(self):
        """Initializes a logger instance"""
        logging.basicConfig(level=logging.DEBUG,
                format='(%(asctime)s)[%(levelname)s] %(message)s',
                filename=os.path.join(self.output_dir, 'eve.log'))

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
        parser.add_argument('-f', '--fasta', required=True,
                            help='Location of genome sequence file to use.')
        parser.add_argument('-g', '--gff', required=True,
                            help='Location of GFF annotation file to use.')
        parser.add_argument('-m', '--mapper', default='bwa',
                            help='Mapper to use for read alignment')
        parser.add_argument('-n', '--num-threads', default='4',
                            help='Maximum number of threads to use')
        parser.add_argument('-t', '--training-set',
                            help='Run EVE in training mode')
        parser.add_argument('-d', '--variant-detectors',
                            default='gatk,mpileup,varscan',
                            help=('Comma-separated list of the variant '
                                  'detectors to be used.'))
        parser.add_argument('-o', '--output-directory',
                            default='output/{timestamp}',
                            help=('Location to store intermediate and output '
                                  'files'))
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

    def load_config(self):
        """Loads the general configuration file for EVE"""
        self.config = configparser.ConfigParser()
        self.config.read('config/eve.cfg')

if __name__ == '__main__':
    app = EVE(sys.argv)
    sys.exit(app.run())
