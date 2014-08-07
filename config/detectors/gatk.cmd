java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R {reference} -I {bam} -o {vcf_unfiltered} -nt {threads}
java -jar GenomeAnalysisTK.jar -T VariantFiltration -R {reference} -V {vcf_unfiltered} -o {vcf_filtered} --filterExpression "DP < 5" --filterName "DepthFilter"
