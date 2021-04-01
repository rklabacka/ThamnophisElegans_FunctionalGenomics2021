# Population Genetic Analysis

In the sub-folders in datafiles are the outputs from the shell scripts and the outputs from the R script (provided in datafiles folder). The steps done in shell are outlined in text below. For the most part, the following github repo was used heavily: https://github.com/simonhmartin/genomics_general.

## Per Gene Analysis:

First, we conducted analysis at the level of gene, including all SNPs that were in the final callset. 

1. Using the python script "parseVCF.py" from the [Genomics General repository](https://github.com/simonhmartin/genomics_general) linked above, the VCF file was first converted to a "geno" formatted file. 

Code used was as follows:

```sh
python $path/VCF_processing/parseVCF.py -i $vcf --skipIndels | gzip > $geno
```

2. The python script "popgenWindows.py" from the [Genomics General repository](https://github.com/simonhmartin/genomics_general) was used to estimate population genetic diversity and divergene in each gene region. FST and DXY were estimated using two different population assignments. The first was at the ecotype level (lake vs. meadow) and the second was at the level of each individual population (N=12). 

Code used was as follows:

```sh
python $path/popgenWindows.py --popsFile $pops --windCoords $bed -g $geno -o genes.windows.csv.gz -f phased -m 1 -T 4 --windType predefined --writeFailedWindows -p MAH -p MER -p PVM -p SUM -p STO -p CHR -p RON -p ROC -p ELF -p NAM -p MAR -p PAP

python $path/popgenWindows.py --popsFile $pops2 --windCoords $bed -g $geno -o genes.windows.meadow-v-lake.csv.gz -f phased -m 1 -T 4 --windType predefined --writeFailedWindows -p Meadow -p Lake
```

3. The gene IDs were added to the output files.

Code used was as follows:

```sh
cat header.txt SeqCap_CapturedGenes.bed >geneIDs.txt
awk '{print $4}' geneIDs.txt >geneIDs.only.txt 
paste genes.windows.meadow-v-lake.csv geneIDs.only.txt -d, > genes.meadow-v-lake.csv
paste genes.windows.csv geneIDs.only.txt -d, > genes.all-pops.csv
```

4. The files from shell were read into R to analyze further (see [R script here](https://github.com/rklabacka/ThamnophisElegans_FunctionalGenomics2021/blob/main/pop_gen_stats/datafiles/pop_gen_stats_FINAL.R)). The main outputs are "fst_within_v_between.csv" and "dxy_within_v_between.csv". This file summarizes the comparison of fst/dxy. First, fst among all pairwise lake vs lake (l_v_l) estimates and meadow vs meadow (m_v_m) estimates was computed by t-test in R. The p-value is reported in the column "within". Similarly, the l_v_l and m_v_m estimates were combined and compared again lake vs meadow (l_v_m) estimates in a t-test in R. The p-value of this test is reported in the column "between". The mean of each group of estimates is reported in the subsequent columns. Additionally, the full FST estimate between all lake samples vs all meadow sample agnostic of population is also provided for comparison. GeneID is on the INFO column. This was repeated for DXY and provided in the second file.


## Per Site Analysis:

Second, we conducted an analysis at the site level, including only SNPs that resulted in a change in the amino acid of the protein sequence (nonsynonymous variants). For the shell code, the biggest difference was in Step #2 where instead of using predefined gene region windows, the --windType argument "sites" was used to estimate a per site diversity and divergence measure for each site. 

Additionally, the population allele count information was obtained and summarized in R. First, using vcftools, the individual population samples were recoded from the original VCF file. This was done using a loop with an array of population IDs.

Code used was as follows:

```sh
for i in "${pops[@]}"
do
    #make list of samples for this population
    awk -v pop=$i '$2==pop {print $1}' pop_list.txt > $i.txt
    wc -l $i.txt

    #split out samples for this pop into a separate VCF file to get population specific allele counts
    vcftools --vcf $vcf --out $i --counts --keep $i.txt
done
```

Later, the tool "vcf-compare" in vcftools was used to get the information used for the upset plot in R. This was done for both the per site and per gene VCF files.

Code used was as follows:

```sh
vcftools --vcf $vcf --recode --out Lake --keep lake.txt --min-alleles 2 --non-ref-ac 1
vcftools --vcf $vcf --recode --out Meadow --keep meadow.txt --min-alleles 2 --non-ref-ac 1

#compare vcfs
bgzip -f Lake.recode.vcf
bgzip -f Meadow.recode.vcf
tabix -f -p vcf Lake.recode.vcf.gz
tabix -f -p vcf Meadow.recode.vcf.gz
vcf-compare Lake.recode.vcf.gz Meadow.recode.vcf.gz >vcf.compare.txt
```
