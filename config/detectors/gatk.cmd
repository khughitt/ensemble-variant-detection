gatk -T UnifiedGenotyper -R {reference} -I {bam} --dbsnp dbsnp.vcf -o {outfile} -nt {threads} -stand_call_conf 30 -stand_emit_conf 10 -dcov 250
