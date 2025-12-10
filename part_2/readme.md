## Pipline in python
In this section, the pipeline is write in a python scripty `pipeline_for_merging_results_from_two_variant_callers.py`. Two variant callers, bcftools and snippy, are used, and the output vcf files from them are merged into one vcf file, for example `Ecoli_simulated_merged.vcf.gz`, and extracted the variants and positions in a csv file. There are three tools applied in the pipeline, minimap2 is used for mapping reads, samtools is used for converting/sorting the output into a bam file, snippy is used for variant calling with a minimal read depth set to 10, and bcftools is used as a variant caller and combining two vcf files. The samtools can also be used for looking at pileup evidence using command line `samtools index Ecoli_simulated_mapped_reads.bam && samtools tview Ecoli_simulated_mapped_reads.bam EcoliK12-MG1655.fasta`.


The input files of the pipeline are, for example using the simulated genome:
    `reference "EcoliK12-MG1655.fasta"`,
    `read_r1 "EcoliK12_simulated_read_R1.fastq"`,
    `read_r2 "EcoliK12_simulated_read_R2.fastq"`,

and the output files and directory are:
    `bam_file "Ecoli_simulated_mapped_reads.bam"`,
    `vcf_bcf "Ecoli_simulated_variants_bcf.vcf"`,
    `vcf_bcf_gz "Ecoli_simulated_variants_bcf.vcf.gz"`,
    `snippy_dir "Ecoli_simulated_snippy_results"`,
    `snippy_vcf_gz "Ecoli_simulated_snippy_results/snps.vcf.gz"`,
    `merged_vcf "Ecoli_simulated_merged.vcf.gz"`.

The result from merged vcf file is extracted into a csv fils based on the GT (genotype) label in the vcf files. The variants can not be detected by the variant caller are labelled as missing. The final output is checked with tview and IGV.

## Checking low confident variants
`count_low_confident_variant.py` could be used for checking the variants with low confident by screening them with low QUAL score, fake heterozygous, or mismatched in two variant callers. 

## Output directory
The reference genomes, real reads, and the output of code in this section can be found in directory: https://jhub.climb.ac.uk/hub/user-redirect/lab/tree/yanyan/Task_2/part_2
