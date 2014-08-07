"""
Variant Detector Classes

TODO: parallelize execution of different variant callers.
"""
import logging
import subprocess

class VariantDetector(object):
    """Base Detector class"""
    def __init__(self, bam, fasta, conf, working_dir, threads):
        """Create a detector instance"""
        self.commands = self.parse_command_template(conf)
        self.bam = bam
        self.fasta = fasta
        self.working_dir = working_dir
        self.threads = threads

    def parse_command_template(self, filepath):
        """Parses a configuration file containing options for the variant
           detector"""
        return [x.strip() for x in open(filepath).readline()]

    def run(self):
        """Runs the given detectors"""
        pass

class GATKDetector(VariantDetector):
    """
    GATKDetector variant detector
    """
    def __init__(self, bam, fasta, conf, working_dir, threads):
        super().__init__(bam, fasta, conf, working_dir, threads)

    def run(self):
        """Run GATK"""
        # GATK otput filepath
        outfile = os.path.join(self.working_dir, 'vcf', 'gatk.vcf')

        cmd = self.commands[0].format(
            reference=self.fasta, bam=self.bam, outfile,
            threads=self.max_threads
        )

        logging.debug(cmd)
        subprocess.call(cmd, shell=True)


class MpileupDetector(VariantDetector):
    """
    Mpileup variant detector

    This class interfaces with the SAMtools Mpileup utility for detecting
    ___ variants.
    """
    def __init__(self, bam, fasta, conf, working_dir, threads):
        super().__init__(bam, fasta, conf, working_dir, threads)

    def run(self):
        """Run the Mpile detection pipeline"""
        # Part 1: mpileup
        bcf_output = os.path.join(self.working_dir, 'vcf', 'var.raw.bcf')
        cmd1 = self.commands[0].format(fasta=self.fasta, bam=self.bam,
                                       bcf_output=bcf_output)
        logging.debug(cmd1)
        subprocess.call(cmd1, shell=True)

        # Part 2: vcfutils
        outfile = os.path.join(self.working_dir, 'vcf', 'mpileup.vcf')
        cmd2 = self.commands[1].format(bcf_output=bcf_output, output=outfile)
        logging.debug(cmd2)
        subprocess.call(cmd2, shell=True)

