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
import pandas
import logging
import argparse
import datetime
import platform
import subprocess
import configparser
import numpy as np
from csv import DictReader
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder,Imputer
from sklearn.externals import joblib
from eve import detectors,mappers

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


class EVE(object):
    """Ensemble Variant Detection"""
    def __init__(self, argv):
        # parse arguments
        self.args = self.parse_args(argv)

        # determine working/output directory to use
        if self.args.output_dir == 'output/{timestamp}':
            now = datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
            self.output_dir = self.args.output_dir.format(timestamp=now)
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

        # split fastq reads into two variables
        (reads1, reads2) = self.args.input_reads

        # output file
        prefix = os.path.basename(
                    os.path.commonprefix([reads1, reads2])).strip("_")
        sam_filepath = "aln_%s.sam" % prefix

        # filepath for mapped reads
        mapped_reads = os.path.join(self.output_dir, 'mapped', sam_filepath)

        # load mapper
        if 'bam' not in self.args:
            if not os.path.exists(mapped_reads):
                self.mapper = mappers.BWAMemMapper(self.args.fasta, reads1, reads2,
                                                mapped_reads, self.args.num_threads)
            else:
                self.args.bam = mapped_reads.replace('.sam', '.bam')

    def run(self):
        """Main application process"""
        # map reads
        if hasattr(self, 'mapper'):
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
            clf_filepath = os.path.join(self.output_dir, 'random_forest.pkl')

            # check to see if classifier already exists
            # @TODO encode parameter information in filepath to ensure same
            # settings were used
            #if not os.path.exists(clf_filepath):
            training_df = self.build_training_set(df)

            # cast fields to string in case results come directly from PyVCF
            # in which case the values would be VCF object instances (e.g.
            # vcf.model._Substitution.
            for x in ['gatk_filtered', 'mpileup', 'varscan_snps']:
                training_df[x] = training_df[x].astype(str)

            # training
            (classifier, training_set, test_set, features, target_classes) = (
                self.train_random_forest(training_df, clf_filepath)
            )
            #else:
            #    classifier = joblib.load(clf_filepath)

        # perform prediction
        self.predict_variants(classifier, test_set, features,
                              target_classes)

        # output final VCF

    def predict_variants(self, classifier, test_set, features, target_classes):
        """Uses a trained Random Forest classifier to predict variants"""
        cls_predict = classifier.predict(test_set[features])

        # map integer predictions back to the original classes
        predictions = target_classes[cls_predict]
        actual = target_classes[np.array(test_set['actual'])]

        pandas.crosstab(actual, predictions,
                        rownames=['actual'], colnames=['preds'])

        # variable importantance
        # http://nbviewer.ipython.org/github/rauanmaemirov/kaggle-titanic101/blob/master/Titanic101.ipynb
        feature_importance = classifier.feature_importances_
        feature_importance = 100.0 * (feature_importance / feature_importance.max())

        sorted_idx = np.argsort(feature_importance)
        pos = np.arange(sorted_idx.shape[0]) + .5

        plt.barh(pos, feature_importance[sorted_idx], align='center')
        plt.yticks(pos, test_set[features].columns[sorted_idx])
        plt.xlabel('Relative feature importance')
        plt.ylabel('Feature name')
        plt.title('EVE Random Forest Variable Importance');
        plt.savefig(os.path.join(self.output_dir,
                    'EVE_Variable_Importance.png'), bbox_inches='tight')

    def train_random_forest(self, df, clf_filepath):
        """
        Trains a random forest classifier on the specified DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            A pandas DataFrame containing both the feature and target columns
            to use for training.
        clf_filepath : str
            Filepath to store classifier at after training.

        Returns
        -------
        classifier : sklearn.ensemble.RandomForestClassifier
            A trained scikit.learn RandomForestClassifier instance instance.
        training_set : DataFrame
            Training set DataFrame
        test_set : DataFrame
            Test set DataFrame
        features : list
            List of feature names used for classification
        target_classes : list
            List of target classes

        TODO
        ----
        * Generalize code to handle arbitrary features.
        * Cross-validation
          (http://scikit-learn.org/stable/modules/cross_validation.html)

        References
        ----------
        |http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
        |http://scikit-learn.org/stable/tutorial/basic/tutorial.html#model-persistence
        |http://blog.yhathq.com/posts/random-forests-in-python.html

        """
        # Encode categorical features using one hot encoding
        features = []

        for orig in ['gatk_filtered', 'mpileup', 'varscan_snps']:
            # replace NaNs to be consistent with target
            df[orig] = df[orig].replace(float('nan'), 'X')

            one_hot = pandas.get_dummies(df[orig])
            for val in one_hot.columns:
                new = "%s=%s" % (orig, val)
                df[new] = one_hot[val]
                features.append(new)

            # delete original column
            df = df.drop(orig, 1)

        # add quality scores and depth to the list of features to test
        features = features + ['depth',  'gatk_filtered_qual',  'mpileup_qual',
                               'varscan_snps_qual']

        # replace nans and encode categorical variable as integers
        df.actual = df.actual.replace(float('nan'), 'X')
        df.actual = df.actual.astype(str)

        encoder = LabelEncoder()
        df.actual = encoder.fit_transform(df.actual)

        #encoder.fit(df.actual)
        #df.actual = encoder.transform(df.actual)

        # impute missing data for quality scores
        for x in ['gatk_filtered_qual',  'mpileup_qual', 'varscan_snps_qual']:
            imputer = Imputer(axis=1)
            df[x] = imputer.fit_transform(df[x])[0]

        # split into training / test data
        df['is_train'] = np.random.uniform(0, 1, len(df)) <= .80
        training_set = df[df['is_train'] == True]
        test_set = df[df['is_train'] == False]

        # train the classifier
        classifier = RandomForestClassifier(n_jobs=-1)
        classifier.fit(training_set[features], training_set['actual'])

        # store classifier
        joblib.dump(classifier, clf_filepath)

        # get the original target classes
        num_classes = len(set(df.actual))

        # array(['A', 'C', 'G', 'T', 'X'], dtype=object
        target_classes = encoder.inverse_transform(range(num_classes))

        return (classifier, training_set, test_set, features, target_classes)

    def build_training_set(self, df):
        """Adds actual values to the end of the combined dataset"""
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier

        # IUPAC codes (needed to parse wgsim output)
        iupac = {
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['G', 'C'],
            'W': ['A', 'T'],
            'K': ['G', 'T'],
            'M': ['A', 'C'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'G', 'C', 'T']
        }

        # Add "truth" values
        df['actual'] = np.repeat(float('nan'), df.shape[0])

        # wgsim
        if (self.args.wgsim):
            fp = open(self.args.training_set)
            wgsim_fields = ['chr', 'pos', 'before', 'after', 'haplotype']
            reader = DictReader(fp, delimiter='\t', fieldnames=wgsim_fields)

            for row in reader:
                # skip insertions, etc.
                if row['before'] not in ['A', 'G', 'C', 'T']:
                    continue
                #if row['after'] not in iupac.keys():
                if row['after'] not in ['A', 'G', 'C', 'T']:
                    continue
                pos = int(row['pos'])
                if pos in df.index:
                    df.loc[df[df.index == pos].index, 'actual'] = row['after']

        else:
            # If training set does not come from wgsim, for now, assume
            # it is a VCF from Genome in a Bottle...
            # load VCF containing true answers
            reader = vcf.Reader(open(self.args.training_set))

            for record in reader:
                #if record.POS in df['position'].values:
                if record.POS in df.index:
                    # there is probably a cleaner way to do this, but I can't
                    # think of it right now...
                    df.loc[df[df.index == record.POS].index, 'actual'] = record.ALT[0]

        df.to_csv(os.path.join(self.output_dir, "combined_training_set.csv"),
                  index_label='position')

        return df

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
        return pandas.DataFrame.from_dict(combined_dict)

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
        parser.add_argument('--wgsim', action='store_true',
                            help='Use wgsim output for training')
        parser.add_argument('-d', '--variant-detectors',
                            default='gatk,mpileup,varscan',
                            help=('Comma-separated list of the variant '
                                  'detectors to be used.'))
        parser.add_argument('-o', '--output-dir',
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
