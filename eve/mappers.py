"""
Read-mapping classes
"""
import os
import logging
import subprocess

class Mapper(object):
    """Base read mapper class"""
    def __init__(self, reference, fastq1, fastq2, outfile, max_threads):
        """Create a detector instance"""
        self.reference = reference
        self.fastq1 = fastq1
        self.fastq2 = fastq2
        self.outfile = outfile
        self.max_threads = max_threads

    def run(self, args):
        """Runs the given mappers"""
        #process = subprocess.Popen(args,
        #                           stdout=subprocess.PIPE,
        #                           stderr=subprocess.PIPE,
        #                           shell=True)
        #stdout, stderr = process.communicate()

        # Going to cheat for a bit...
        subprocess.call(" ".join(args), shell=True)

        #if stdout:
        #    logging.info(stdout)
        #if stderr:
        #    logging.error(stderr)

        #return process.returncode
        return self.outfile

class BWAMemMapper(Mapper):
    """Burrows-Wheeler Aligner Mapper class"""
    def __init__(self, reference, fastq1, fastq2, outfile, max_threads):
        super().__init__(reference, fastq1, fastq2, outfile, max_threads)

    def run(self):
        """Run BWA mapping command"""
        cmd1 = "bwa mem -t {threads} {reference} {fastq1} {fastq2} > {output}".format(
                    reference=self.reference,
                    fastq1=self.fastq1, fastq2=self.fastq2,
                    threads=self.max_threads, output=self.outfile
        )
        logging.debug(cmd1)

        # Call base Detector.run method
        #args_list = cmd.split(' ')
        #super().run(args_list)
        subprocess.call(cmd1, shell=True)

        # Convert to BAM
        bam = self.outfile.replace(".sam", ".bam")

        cmd2 = "samtools view -bS {sam} > {bam}".format(
            sam=self.outfile, bam=bam
        )

        logging.debug(cmd2)
        subprocess.call(cmd2, shell=True)

        # Sort and index
        bam_sorted = bam.replace('.bam', '_sorted')
        cmd3 = "samtools sort {bam} {bam_sorted}".format(
            bam=bam, bam_sorted=bam_sorted
        )
        logging.debug(cmd3)
        subprocess.call(cmd3, shell=True)

        # Add read groups
        bam_sorted_rg = "%s_RG.bam" % bam_sorted
        cmd4 = ("java -jar AddOrReplaceReadGroups.jar I={bam_sorted}.bam "
                "O={bam_sorted_rg} RGLB=lib RGPL=illumina "
                "RGPU=4410 RGSM=Project").format(
            bam_sorted=bam_sorted, bam_sorted_rg=bam_sorted_rg
        )
        logging.debug(cmd4)
        subprocess.call(cmd4, shell=True)

        # Index BAM file
        cmd5 = "samtools index {bam_sorted_rg}".format(
            bam_sorted_rg=bam_sorted_rg
        )

        logging.debug(cmd5)
        subprocess.call(cmd5, shell=True)

        return bam_sorted_rg

