"""
Variant Detector Classes

TODO: parallelize execution of different variant callers.
"""
import logging
import subprocess
from yaml import load, Loader

class VariantDetector(object):
    """Base Detector class"""
    def __init__(self, bam, fasta, conf, working_dir):
        """Create a detector instance"""
        self.config = self.parse_config(conf)
        self.bam = bam
        self.fasta = fasta
        self.working_dir = working_dir

    def parse_config(self, filepath):
        """Parses a YAML configuration file containing options for the variant
           detector"""
        if filepath.endswith('.yaml'):
            return load(open(filepath), Loader=Loader)
        elif filepath.endswith('.txt'):
            return [x.strip() for x in open(filepath).readline()]

    def get_arg_list(self):
        """Returns a list of command-line arguments for the detector"""
        args = [self.config['command']]

        for key in self.config['parameters']:
            value = self.config['parameters'][key]
            if value is None:
                args.append(" %s" % key)
            elif key.startswith('--'):
                args.append((" %s=%s" % (key, value)))
            else:
                args.append((" %s %s" % (key, value)))

        return(args)

    def run(self):
        """Runs the given detectors"""
        logging.debug(" ".join(self.get_arg_list()))
        print(" ".join(self.get_arg_list()))

        process = subprocess.Popen(self.get_arg_list(),
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)

        stdout, stderr = process.communicate()

        if stdout:
            logging.info(stdout)
        if stderr:
            logging.error(stderr)

        return process.returncode

class GATKDetector(VariantDetector):
    """
    GATKDetector variant detector
    """
    def __init__(self, bam, fasta, conf, working_dir):
        super().__init__(bam, fasta, conf, working_dir)

class MpileupDetector(VariantDetector):
    """
    Mpileup variant detector

    This class interfaces with the SAMtools Mpileup utility for detecting
    ___ variants.
    """
    def __init__(self, bam, fasta, conf, working_dir):
        super().__init__(bam, fasta, conf, working_dir)

    def run(self):
        """Run the Mpile detection pipeline"""
        # Part1: mpileup
        bcf_output = os.path.join(self.working_dir, 'vcf', 'var.raw.bcf')
        cmd1 = self.config[0].format(fasta=self.fasta, bam=self.bam,
                                     bcf_output=bcf_output)
        print(cmd1)

        # hard-coded for now...
        # samtools mpileup -uf /scratch/summerschool/Data/Reference_Genomes/hsa_37.p5/hs_ref_GRCh37.p5_chr22.fa
        # /scratch/summerschool/Data/Data_From_Tuebingen/NA12878_01.bam \
        # | bcftools view -bvcg - \
        # > var.raw.bcf \


