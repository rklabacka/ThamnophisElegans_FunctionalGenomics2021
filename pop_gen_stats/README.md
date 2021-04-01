# Population Genetic Analysis

## Per Gene Analysis:

First, we conducted analysis at the level of gene, including all SNPs that were in the final callset. 

FST and DXY were estimated using two different population assignments. The first was at the ecotype level (lake vs. meadow) and the second was at the level of each individual population (N=12). For the all populations, a comparison of fst among all lake vs lake (l_v_l) estimates and meadow vs meadow (m_v_m) estimates was computed by t-test in R. The p-value is reported in the column "within". Similarly, the l_v_l and m_v_m estimates were combined and compared again lake vs meadow (l_v_m) estimates in a t-test in R. The p-value of this test is reported in the column "between". The mean of each group of estimates is reported in the subsequent columns. Additionally, the full FST estimate between all lake samples vs all meadow sample agnostic of population is also provided for comparison. GeneID is on the INFO column.

In the subsequent folder are the outputs from the shell scripts and the outputs from the R script. All scripts are included.


## Per Site Analysis:

Second, we conducted an analysis at the site level, including only SNPs that resulted in a change in the amino acid of the protein sequence (nonsynonymous variants).
