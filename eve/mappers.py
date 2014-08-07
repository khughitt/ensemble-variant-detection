"""
Read-mapping classes
"""
import os
import logging
import subprocess

class Mapper(object):
    """Base read mapper class"""
    def __init__(self, fasta1, fasta2, outfile):
        """Create a detector instance"""
        self.fasta1 = fasta1
        self.fasta2 = fasta2
        self.outfile = outfile

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
    def __init__(self, fasta1, fasta2, outfile):
        super().__init__(fasta1, fasta2, outfile)

    def run(self):
        """Run BWA mapping command"""
        # @TODO: accept number of threads as argument
        cmd = "bwa mem -t {threads} {fasta1} {fasta2} > {output}".format(
                    fasta1=self.fasta1, fasta2=self.fasta2,
                    threads=32, output=self.outfile
        )
        logging.debug(cmd)

        # Run BWA
        super().run(cmd.split(' '))

