samtools mpileup -f {reference} {bam} > {mpileup_output}
java -jar VarScan.v2.3.7.jar mpileup2snp {mpileup_output} --min-coverage 5 --output-vcf 1 > {varscan_snps}
java -jar VarScan.v2.3.7.jar mpileup2indel {mpileup_output} --min-coverage 5 --output-vcf 1 > {varscan_indels}
