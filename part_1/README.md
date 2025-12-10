## Step_1: preparing simulated mutated genome
In this step, 300 SNPs and 20 Indels (1 - 10 bases) are generated using script `make_mutation_genome.py`. 
The original input, i.e., reference genome sequence, must be provided as a fasta file. It will be converted into a sequence used to construct the mutated genome.


For the code, first of all, 300 SNPs 
are randomly picked and replace the bases with 4 special characters. The insertions are also protected by wrapping up with insertion markers. In this case, while doing the deletion, it will avoid any special characters and wrapped ranges, and will not change the positions of either insertions or the SNPs. After all operations indels, the special characters of picked SNPs are replaced by random choice of 3 different bases, and all wrapping markers are removed to get the final simulated mutated genome.
The simulated mutated genome is saved in`simulated_mutate_genome.txt` for later usage. All simulated mutations are recorded in a CSV file named `simulated_mutation.csv` with four columns: Operation (SNP, Insertion, Deletion), POS (position), REF (base in reference genome), and ALT (simulated mutation).

## Step_2: simulating paired-end reading from Illumina and stored outcome in two fastq files
In this step, the `simulated_mutated_genome.txt ` file generated in step_1 is used to prepare 100-bases paired-end reads with a 30 reading depth from Illumina. The resulting reads are saved as a pair of fastq files: `simulated_read_1.fastq` and `simulated_read_2.fastq`. These fastq files will be used to map the reads to the reference genome and to run bcftools for variant calling. The mapping and variant calling processes are performed on the server using command line.

## Step_3: comparing the result from variant caller with simulated mutations record
The outcome of variant calling, `vcf_variants.vcf`, is downloaded from the server and compared with the recorded `simulated_mutation.csv`. The output, `merged_result.csv`, contains six columns, POS (position), REF (base in reference genome), ALT (base in mutated genome), Type (SNP or INDEL), Match_Status (MATCH or MISMATCH), adn Source (Both, CSV_only, and VCF_only).

## Output directory
The reference genomes used and the output, of the code in this section can be found in directory: https://jhub.climb.ac.uk/hub/user-redirect/lab/tree/yanyan/Task_2/part_1.
