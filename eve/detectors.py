"""
Variant Detector Classes

TODO: parallelize execution of different variant callers.
"""
import os
import logging
import subprocess

class VariantDetector(object):
    """Base Detector class"""
    def __init__(self, bam, fasta, conf, output_dir, threads, location):
        """Create a detector instance"""
        self.commands = self.parse_command_template(conf)
        self.bam = bam
        self.fasta = fasta
        self.output_dir = output_dir
        self.threads = threads
        self.location = location

    def parse_command_template(self, filepath):
        """Parses a configuration file containing options for the variant
           detector"""
        with open(filepath) as fp:
            return [x.strip() for x in fp.readlines()]

    def run(self):
        """Runs the given detectors"""
        pass

class GATKDetector(VariantDetector):
    """
    GATKDetector variant detector
    """
    def __init__(self, bam, fasta, conf, output_dir, threads, location):
        super().__init__(bam, fasta, conf, output_dir, threads, location)

    def run(self):
        """Run GATK"""
        logging.info("Running GATK")

        # GATK output filepaths
        unfiltered_vcf = os.path.join(self.output_dir, 'vcf', 'gatk_unfiltered.vcf')
        filtered_vcf = os.path.join(self.output_dir, 'vcf', 'gatk_filtered.vcf')

        # If output files already exist, stop here
        if os.path.exists(filtered_vcf):
            logging.info("GATK output already exists. Skipping...")
            return filtered_vcf

        # Find all SNPs and indels, regardless of coverage
        cmd = self.commands[0].format(
            jar=self.location, reference=self.fasta, bam=self.bam,
            vcf_unfiltered=unfiltered_vcf, threads=self.threads
        )

        logging.debug(cmd)
        subprocess.call(cmd, shell=True)

        # Filter out low-coverage hits
        cmd = self.commands[1].format(
            jar=self.location, reference=self.fasta,
            vcf_unfiltered=unfiltered_vcf, vcf_filtered=filtered_vcf
        )

        logging.debug(cmd)
        subprocess.call(cmd, shell=True)

        # Clean up
        os.unlink(unfiltered_vcf)
        os.unlink(unfiltered_vcf + ".idx")

        return filtered_vcf

class MpileupDetector(VariantDetector):
    """
    Mpileup variant detector

    This class interfaces with the SAMtools Mpileup utility for detecting
    ___ variants.
    """
    def __init__(self, bam, fasta, conf, output_dir, threads, location):
        super().__init__(bam, fasta, conf, output_dir, threads, location)

    def run(self):
        """Run the Mpile detection pipeline"""
        # Part 1: mpileup
        logging.info("Running Mpileup")

        # output filepaths
        bcf_output = os.path.join(self.output_dir, 'vcf', 'var.raw.bcf')
        vcf_output = os.path.join(self.output_dir, 'vcf', 'mpileup.vcf')

        # If output files already exist, stop here
        if os.path.exists(vcf_output):
            logging.info("Mpileup output already exists. Skipping...")
            return vcf_output

        # Otherwise run Mpileup
        cmd1 = self.commands[0].format(fasta=self.fasta, bam=self.bam,
                                       bcf_output=bcf_output)
        logging.debug(cmd1)
        subprocess.call(cmd1, shell=True)

        # Part 2: vcfutils
        cmd2 = self.commands[1].format(bcf_output=bcf_output,
                                       output=vcf_output)
        logging.debug(cmd2)
        subprocess.call(cmd2, shell=True)

        # Cleanup
        os.unlink(bcf_output)

        return vcf_output

class VarScanDetector(VariantDetector):
    """
    VarScanDetector variant detector
    """
    def __init__(self, bam, fasta, conf, output_dir, threads, location):
        super().__init__(bam, fasta, conf, output_dir, threads, location)

    def run(self):
        """Run VarScan"""
        logging.info("Running VarScan")

        # VarScan intermediate and output filepaths
        mpileup_outfile = os.path.join(self.output_dir, 'mpileup',
                                       'varscan.mpileup')
        varscan_snps_vcf = os.path.join(self.output_dir, 'vcf',
                                       'varscan_snps.vcf')
        varscan_indels_vcf = os.path.join(self.output_dir, 'vcf',
                                          'varscan_indels.vcf')

        # If output files already exist, stop here
        if os.path.exists(varscan_snps_vcf):
            logging.info("VarScan output already exists. Skipping...")
            return varscan_snps_vcf

        # Step 1: mpileup
        cmd1 = self.commands[0].format(
            reference=self.fasta, bam=self.bam, mpileup_output=mpileup_outfile
        )

        logging.debug(cmd1)
        subprocess.call(cmd1, shell=True)

        # Step 2: pileup2snp
        cmd2 = self.commands[1].format(
            jar=self.location, mpileup_output=mpileup_outfile,
            varscan_snps=varscan_snps_vcf
        )

        logging.debug(cmd2)
        subprocess.call(cmd2, shell=True)

        # Step 3: pileup2indel
        # (skipping this for now to speed things up...)
        #cmd3 = self.commands[2].format(
        #    jar=self.location, mpileup_output=mpileup_outfile,
        #    varscan_indels=varscan_indels_vcf
        #)

        #logging.debug(cmd3)
        #subprocess.call(cmd3, shell=True)

        #return [varscan_snps_vcf, varscan_indels_vcf]
        return varscan_snps_vcf

