## Discussion: What is the precision and recall of the SNP caller?
There are 320 mutations in the simulated genome and bcftools can detect 319 of them. In the `merged_result.csv`, two types of mismatch can be observed.

![Simulated SNP cannot be detected 
by bcftools](image/SNP_CSV_only.png)
**Figure 1.** Example of simulated SNP cannot be detected by bcftools, but recorded in the csv file.


![Indel can be detected by bcftools, 
but record in a different way with 
cvs file](image/Indel_VCF_and_CSV.png)
**Figure 2.** Example of indel can be detected by bcftools, but recorded in a different way with cvs file.

In the simulated mutated genome, there is only one SNP that cannot be detected by bcftools (Figure 1), and all other mismatches of indels are due to differences in recording methods. Based on the example in Figure 2, the simulated mutation is located at position 861277, with the recorded reference base as "G" at position 861276. The result for bcftools also starts from position 861276, but the recorded reference is "GT". This discrepancy can be explained by the underlying algorithm of bcftools. This mutation can be recorded in different ways by the variant caller; bcftools follows the parsimony principle and records the mutation in a different way rather than simply including several bases insertion at a certain position.

In the test run of _E. coli K12_, bcftools detected 299 SNPs and 20 indels. In 20 indels, 7 of them are labelled as mismatch but because of the different way of recording as described before. In this case, for _E. coli K12_, the rate of recall is $（319/320）* 100\% = 99.69\%$, and precision is same.

In the test run of _P. falciparum_, bcftools detected 300 SNPs with one false detection and one missing detection, and 20 indels, where 12 of them labeled as mismatch but can be explained same as above. Therefore, the rate of recall is $（318/320） * 100\% = 99.38\%$, and precision is same.
