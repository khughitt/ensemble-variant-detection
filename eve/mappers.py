"""
Read-mapping classes
"""
class Mapper(object):
    """Base read mapper class"""
    def __init__(self, fasta1, fasta2, conf, working_dir, outfile_mapping):
        """Create a detector instance"""
        self.config = self.parse_config(conf)
        self.fasta1 = fasta1 
        self.fasta2 = fasta2
        self.working_dir = working_dir
        self.outfile_mapping = outfile_mapping            

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
    def __init__(self, fasta1, fasta2, conf, working_dir, outfile_mapping):
        super().__init__(fasta1, fasta2, conf, working_dir, outfile_mapping)
        
    def run(self):
        
        bwa_output = os.path.join(self.working_dir, 'aln-pe.sam')
        cmd = "bwa mem -t {threads} {fasta1} {fasta2} > {output}".format(fasta1=self.fasta1, 
                        fasta2=self.fasta2, threads=self.threads, output=bwa_output)      
        print(cmd)
        super().run(cmd.split(' '))                
        # bwa mem -t 32 ../Reference_Genomes/hg19/hg19_complete.fasta NA12878_01_chr22_L001_R1_001.fastq NA12878_01_chr22_L001_R2_001.fastq  > aln-NA12878_01_chr22_L001_untrimmed-pe.sam
        
