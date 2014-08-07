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
        process = subprocess.Popen(args,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if stdout:
            logging.info(stdout)
        if stderr:
            logging.error(stderr)

        return process.returncode

class BWAMemMapper(Mapper):
    """Burrows-Wheeler Aligner Mapper class"""
    def __init__(self, reference, fastq1, fastq2, outfile, max_threads):
        super().__init__(reference, fastq1, fastq2, outfile, max_threads)

    def run(self):
        """Run BWA mapping command"""
        # @TODO: accept number of threads as argument
        cmd = "bwa mem -t {threads} {reference} {fastq1} {fastq2} > {output}".format(
                    reference=self.reference,
                    fastq1=self.fastq1, fastq2=self.fastq2,
                    threads=self.max_threads, output=self.outfile
        )
        logging.debug(cmd)

        # Run BWA
        super().run(cmd.split(' '))

