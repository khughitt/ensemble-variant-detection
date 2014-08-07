samtools mpileup -uf {fasta} {bam} | bcftools view -bvcg - > {bcf_output}
bcftools view {bcf_output} | vcfutils.pl varFilter -d5 -D100 > {output}
