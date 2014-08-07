samtools mpileup -f {reference} {bam} > {mpileup_output}
varscan2 pileup2snp {mpileup_output} > {varscan_snps}
varscan2 pileup2indel {mpileup_output} > {varscan_indels}
