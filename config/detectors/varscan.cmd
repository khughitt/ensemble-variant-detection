samtools mpileup -f {reference} {bam} > {mpileup_output}
java -jar {jar} mpileup2snp {mpileup_output} --min-coverage 5 --output-vcf 1 > {varscan_snps}
java -jar {jar} mpileup2indel {mpileup_output} --min-coverage 5 --output-vcf 1 > {varscan_indels}
