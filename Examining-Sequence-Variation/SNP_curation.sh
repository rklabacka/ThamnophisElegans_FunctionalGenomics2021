#!/bin/sh

function combineDatasets {
cd $WorkingDirectory/variantFiltration
# Jessica's WGS variants at variable sites discovered from q.FullScript_May2020.sh WGS_Genes.recode.vcf.gz WGS_Exons.recode.vcf.gz and WGS_CDS.recode.vcf.gz were copied here frome box already
gunzip SeqCap_Genes.vcf.gz
gunzip SeqCap_Exons.vcf.gz
gunzip SeqCap_CDS.vcf.gz
grep -v "#" SeqCap_Genes.vcf | cut -f 1,2 > SeqCap_Genes.pos 
grep -v "#" SeqCap_Exons.vcf | cut -f 1,2 > SeqCap_Exons.pos 
grep -v "#" SeqCap_CDS.vcf | cut -f 1,2 > SeqCap_CDS.pos 
vcftools --vcf Genes_WGS.recode.vcf --positions SeqCap_Genes.pos --recode --recode-INFO-all --out WGS_Genes
mv WGS_Genes.recode.vcf WGS_Genes.vcf
vcftools --vcf Genes_WGS.recode.vcf --positions SeqCap_Exons.pos --recode --recode-INFO-all --out WGS_Exons 
mv WGS_Exons.recode.vcf WGS_Exons.vcf
vcftools --vcf Genes_WGS.recode.vcf --positions SeqCap_CDS.pos --recode --recode-INFO-all --out WGS_CDS
mv WGS_CDS.recode.vcf WGS_CDS.vcf
bgzip SeqCap_Genes.vcf
bgzip SeqCap_Exons.vcf
bgzip SeqCap_CDS.vcf
bcftools index -f SeqCap_Genes.vcf.gz
bcftools index -f SeqCap_Exons.vcf.gz
bcftools index -f SeqCap_CDS.vcf.gz
bgzip -f WGS_Genes.vcf
bgzip -f WGS_Exons.vcf
bgzip -f WGS_CDS.vcf
bcftools index -f WGS_Genes.vcf.gz
bcftools index -f WGS_Exons.vcf.gz
bcftools index -f WGS_CDS.vcf.gz
echo "merging WGS data to RNA-Seq+Seq-Cap" >> Log.txt
bcftools merge SeqCap_CDS.vcf.gz WGS_CDS.vcf.gz -O v -o Full_CDS.vcf
bcftools merge SeqCap_Exons.vcf.gz WGS_Exons.vcf.gz -O v -o Full_Exons.vcf
bcftools merge SeqCap_Genes.vcf.gz WGS_Genes.vcf.gz -O v -o Full_Genes.vcf
bgzip -f Full_CDS.vcf
bgzip -f Full_Exons.vcf
bgzip -f Full_Genes.vcf
bcftools index -f Full_CDS.vcf.gz
bcftools index -f Full_Exons.vcf.gz
bcftools index -f Full_Genes.vcf.gz
}

function sortVariants {
cd $WorkingDirectory/variantFiltration
cp "$1".vcf.gz "$1"_original.vcf.gz
bcftools index -f "$1".vcf.gz
bcftools norm -d snps "$1".vcf.gz -O v -o "$1"_dupsIDed.vcf
bcftools view -m2 -M2 -v snps "$1"_dupsIDed.vcf > "$1"_dupsRemoved.vcf
bcftools sort "$1"_dupsRemoved.vcf -O z -o "$1"_sorted.vcf.gz
gunzip "$1".vcf.gz
echo "original:  $(grep -v "^#" "$1".vcf | wc -l)" >> Log.txt
echo "with dups removed: $(grep -v "^#" "$1"_dupsRemoved.vcf | wc -l)" >> Log.txt
# rm "$1"_dupsRemoved.vcf
# rm "$1".vcf
mv "$1"_sorted.vcf.gz "$1".vcf.gz
bcftools index -f "$1".vcf.gz
rm "$1".vcf
}

function createPopFiles {
  cd $WorkingDirectory/SNP_analysis
  mkdir -p Populations
  cd Populations
  for sample in `bcftools query -l $WorkingDirectory/variantFiltration/Full_CDS.vcf.gz`
  do
    echo "$sample" >> Samples
  done
  python $pythonScripts/parsePopulations.py Samples
  rm *PIK*.txt
  mkdir -p pairwisePops
  mkdir -p allPops
  mkdir -p pixyPops
  mv *_*.txt pairwisePops
  mv *.pix pixyPops
  cd pixyPops
  rm *PIK*.pix
  cd ..
  mv *.txt allPops
  cd allPops
  rm PIK.txt
  ls *.txt > PopList
  python $pythonScripts/joinPopulations.py PopList ../samples-by-population.txt
} 

function createPairwiseVCFs {
  cd $WorkingDirectory/variantFiltration
  gunzip Full_CDS_missense.vcf.gz
  gunzip Full_CDS_synonymous.vcf.gz
  gunzip Full_Exons.vcf.gz
  cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops
  ls *.txt | cut -d "." -f "1" | sort > PairwisePopsList
  echo -e "PairwiseComparison\tN\tMissense\tSynonymous" >> PairwisePopSegregatingSites.txt
  while read i
  do
    echo "Population $i"
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops
    mkdir -p "$i"/missense "$i"/synonymous "$i"/exons
    mv "$i".txt "$i"/
    bcftools view --samples-file "$i"/"$i".txt --min-ac=1 --no-update \
        $WorkingDirectory/variantFiltration/Full_CDS_missense.vcf > "$i"/missense/Full_CDS_missense_"$i".vcf
    bcftools view --samples-file "$i"/"$i".txt --min-ac=1 --no-update \
        $WorkingDirectory/variantFiltration/Full_CDS_synonymous.vcf > "$i"/synonymous/Full_CDS_synonymous_"$i".vcf
    bcftools view --samples-file "$i"/"$i".txt --min-ac=1 --no-update \
        $WorkingDirectory/variantFiltration/Full_Exons.vcf > "$i"/exons/Full_Exons_"$i".vcf
    n="$(wc -l < "$i")"
    misSNPcount="$(grep -v "^#" "$i"/missense/Full_CDS_missense_"$i".vcf | wc -l)"
    synSNPcount="$(grep -v "^#" "$i"/synonymous/Full_CDS_synonymous_"$i".vcf | wc -l)"
    echo -e "$i\t$n\t$misSNPcount\t$synSNPcount\t$tajD" >> PairwisePopSegregatingSites.txt
    j="$i"
    k="$i"
    cd "$i"/missense/
    bgzip Full_CDS_missense_"$i".vcf
    bcftools index -f Full_CDS_missense_"$i".vcf.gz
    getVCFbyGene CDS Full_CDS_missense_"$i".vcf.gz "$i"
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$j"/synonymous
    bgzip Full_CDS_synonymous_"$j".vcf
    bcftools index -f Full_CDS_synonymous_"$j".vcf.gz
    getVCFbyGene CDS Full_CDS_synonymous_"$j".vcf.gz "$j"
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$k"/exons
    bgzip Full_Exons_"$k".vcf
    bcftools index -f Full_Exons_"$k".vcf.gz
    getVCFbyGene Exons Full_Exons_"$k".vcf.gz "$k"
  done<PairwisePopsList
  cd $WorkingDirectory/variantFiltration
  bgzip Full_CDS_missense.vcf
  bgzip Full_CDS_synonymous.vcf
  bgzip Full_Exons.vcf
}

function getVariantBED {
cd $WorkingDirectory/variantFiltration
gunzip "$1".vcf.gz
vcf2bed < "$1".vcf > $WorkingDirectory/References/"$1"_variants.bed
bgzip "$1"_variants.bed
bgzip "$1".vcf
bcftools index -f bgzip "$1".vcf
}

function functionalAnnotation {
# Step 1: Download and install
  cd ~
  wget wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
  unzip snpEff_latest_core.zip
# Step 2: Create genome annotation database
  cd snpEff
  mkdir -p data
  cd data
  mkdir -p rThaEle1 genomes
  cp $WorkingDirectory/References/SeqCap_CapturedCDS.gff ./rThaEle1/genes.gff
  cp $WorkingDirectory/References/TelagGenome.fasta ./genomes/rThaEle1.fa
  cd ~/snpEff
  # Before running the build command below, make sure you have the configuration file with your genome specified
  # You will also need to make sure you have enough memory alloted for the build
  java -jar /tools/snpeff-5.0/snpEff.jar build -gff3 -v rThaEle1
# Step 3: Run snpEff
  cd $WorkingDirectory/variantFiltration
  gunzip "$1"_CDS.vcf.gz
  java -jar /tools/snpeff-5.0/snpEff.jar -c ~/snpEff/snpEff.config -v rThaEle1 "$1"_CDS.vcf > "$1"_CDS_ann.vcf
  # This didn't work with the raw gff downloaded from genbank (TelagGenome.gff). 
  # Instead, I had to use the "$1"_CapturedCDS.gff file I modified to only include CDS in genes of interest.
# Step 4: Pull out missense and synonymous mutations
  awk '/^#|missense_variant/' "$1"_CDS_ann.vcf > "$1"_CDS_missense.vcf
  awk '/^#|synonymous_variant/' "$1"_CDS_ann.vcf > "$1"_CDS_synonymous.vcf
  echo "total CDS SNP count: $(grep -v "^#" "$1"_CDS_ann.vcf | wc -l)" >> Log.txt
  echo "missense SNP count: $(grep -v "^#" "$1"_CDS_missense.vcf | wc -l)" >> Log.txt
  echo "synonymous SNP count: $(grep -v "^#" "$1"_CDS_synonymous.vcf | wc -l)" >> Log.txt
  bgzip "$1"_CDS.vcf
  bcftools index -f "$1"_CDS.vcf
  bgzip "$1"_CDS_missense.vcf
  bgzip "$1"_CDS_synonymous.vcf
  bcftools index -f "$1"_CDS_missense.vcf.gz
  bcftools index -f "$1"_CDS_synonymous.vcf.gz
  bgzip "$1"_CDS_ann.vcf
  bcftools index -f "$1"_CDS_ann.vcf.gz

}

function getGeneVariants {
# Get CDS SNPs and prepare for extraction
cd $WorkingDirectory/variantFiltration
bgzip -f Full_"$1$2".vcf
bcftools index -f Full_"$1$2".vcf.gz
mkdir -p $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
cp $WorkingDirectory/variantFiltration/Full_"$1""$2".vcf.gz .
bcftools index -f Full_"$1$2".vcf.gz
#+ COMPLETED # Create bed file for each gene
#+ COMPLETED getBEDbyGene $1 $2
# Extract SNPs by gene from vcf
echo "begin SNP extraction by gene from vcf" >> Log.txt
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
getVCFbyGene $1 $WorkingDirectory/variantFiltration/Full_"$1""$2".vcf.gz $2 
}

function getBEDbyGene {
# Create bed file for each gene
cd $WorkingDirectory/References
mkdir -p GeneBEDs
cd GeneBEDs
python $pythonScripts/parseBED.py ../SeqCap_CapturedGenes.bed Full_"$1""$2"_Captures.txt "$1"
sort -u Full_"$1""$2"_Captures.txt > Full_"$1""$2"_CapturedGeneList.txt
}

function getVCFbyGene {
# Extract SNPs by gene from vcf
locivar=0
while read i
  do
  mkdir -p "$i"
  bcftools view -R $WorkingDirectory/References/GeneBEDs/"$i"_"$1".bed "$2" > "$i"/"$i"_"$1$3".vcf
  locusvar="$(grep -v "^#" "$i"/"$i"_"$1$3".vcf| wc -l)"
  echo "$i	$locusvar" >> Full_"$1""$3"_Nvariants.txt
  locivar=$((locusvar + locivar))
done<$WorkingDirectory/References/GeneBEDs/Full_CDS_CapturedGeneList.txt
  echo "total variants: $locivar" >> Full_"$1$3"_Nvariants.txt
}

function getGeneTajD {
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1$2"
echo "Gene  TajD" >> Full_"$1$2"_TajD.txt
while read i
do
  cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1$2"/$i
  variants=($(ls *vcf))
  vcftools --vcf $variants --TajimaD 1000000 --out "$i"
  rm "$i".log
  TajD=$(awk 'FNR == 2 {print $4}' "$i".Tajima.D) 
  TajD=${TajD:="NA"}
  cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1$2"
  echo "$i  $TajD" >> Full_"$1$2"_TajD.txt
done<$WorkingDirectory/References/GeneBEDs/Full_CDS_CapturedGeneList.txt
}

function getPairwisePopGen {
cd $WorkingDirectory
mkdir -p $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats
cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops
echo "begin pairwise pop gen" >> Log.txt
while read i
do
  cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats
  # Add header to gene outfile
  echo -e "Population\tN\tMissenseSNPs\tSynonymousSNPs\tTranscriptSNPs\tTajD" >> "$i".txt
  while read j
  do
      echo "pairwise pops: ${j}" >> $WorkingDirectory/popGenLog.txt
      echo "gene: ${i}" >> $WorkingDirectory/popGenLog.txt
    pop1=`echo $j | awk -F"_" '{print $1}'`
      echo -e "\tPopulation 1: ${pop1}" >> $WorkingDirectory/popGenLog.txt
    pop2=`echo $j | awk -F"_" '{print $2}'`
      echo -e "\tPopulation 2: ${pop2}" >> $WorkingDirectory/popGenLog.txt
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$j"
    n="$(wc -l < "$j".txt)"
      echo -e "\tnumber of individuals: ${n}" >> $WorkingDirectory/popGenLog.txt
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$j"/missense/"$i"
    misSNPcount="$(grep -v "^#" "$i"_CDS_"$j".vcf | wc -l)"
      echo -e "\tmisSNPcount: ${misSNPcount}" >> $WorkingDirectory/popGenLog.txt
    # Use vcftools to calculate Fst for each missense site
    vcftools --vcf "$i"_CDS_"$j".vcf \
        --weir-fst-pop $WorkingDirectory/SNP_analysis/Populations/allPops/"$pop1".txt \
        --weir-fst-pop $WorkingDirectory/SNP_analysis/Populations/allPops/"$pop2".txt \
        --out "$j"_"$i"
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$j"/synonymous/"$i"
    synSNPcount="$(grep -v "^#" "$i"_CDS_"$j".vcf | wc -l)"
      echo -e "\tsynSNPcount: ${synSNPcount}" >> $WorkingDirectory/popGenLog.txt
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$j"/exons/"$i"
    transcriptSNPcount="$(grep -v "^#" "$i"_Exons_"$j".vcf | wc -l)"
      echo -e "\ttranscriptSNPcount: ${transcriptSNPcount}" >> $WorkingDirectory/popGenLog.txt
    # Use vcftools to calculate tajima's D for the exonic regions of each gene
    vcftools --vcf "$i"_Exons_"$j".vcf --TajimaD 1000000 --out "$j"_"$i"
    tajD="$(awk 'NR == 2 {print $4}' $j'_'$i.Tajima.D)"
    tajD=${tajD:="NA"}
      echo -e "\ttajD: ${tajD}" >> $WorkingDirectory/popGenLog.txt
    # Now using genomics-general package for dxy, fst, and pi ---------------------------------
    # Remove this comment and modify header for outfile if you decide to re-implement
    # pixy function to calculate average fst and dxy for exonic regions of each gene
    # pixyPopGen $i $j  
    # cd $WorkingDirectory/SNP_analysis/Populations/pixyPops/$j/$i/
    # fst="$(awk 'NR == 2 {print $6}' "$i"_"$j"_perGene_fst.txt)"
    # fst=${fst:="NA"}
    #   echo -e "\tfst: ${fst}" >> $WorkingDirectory/popGenLog.txt
    # dxy="$(awk 'NR == 2 {print $6}' "$i"_"$j"_perGene_dxy.txt)"
    # dxy=${dxy:="NA"}
    #   echo -e "\tdxy: ${dxy}" >> $WorkingDirectory/popGenLog.txt
    # # the pi values calculated from pixy are average pairwise differences within each population
    # pop1pi="$(awk 'NR == 2 {print $5}' "$i"_"$j"_perGene_pi.txt)"
    # pop1pi=${pop1pi:="NA"}
    #   echo -e "\tpop1pi: ${pop1pi}" >> $WorkingDirectory/popGenLog.txt
    # pop2pi="$(awk 'NR == 3 {print $5}' "$i"_"$j"_perGene_pi.txt)"
    # pop2pi=${pop2pi:="NA"}
    #   echo -e "\tpop2pi: ${pop2pi}" >> $WorkingDirectory/popGenLog.txt
    # --------------------------------------------------------------------------
    # Add values to gene outfile
      echo -e "$j\t$n\t$misSNPcount\t$synSNPcount\t$transcriptSNPcount\t$tajD" >> $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats/"$i".txt
  done<$WorkingDirectory/SNP_analysis/Populations/pairwisePops/PairwisePopsList
  cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats
  python $pythonScripts/addEcotypes.py "$i".txt "$i"_withEcotypes.txt
done<$WorkingDirectory/References/GeneBEDs/Full_CDS_CapturedGeneList.txt
}

function pairwisePopGen2 {
# This function is used to get population genetics statistics from pairwise population comparisons
# It is less expensive (computationally), since the genomics_general package performs the population pairing and statistical calculations
# Fst, Dxy, and Pi are obtained for both windows and sites

# Step 1: Create working environment 
mkdir -p $WorkingDirectory/SNP_analysis/genomics_general
cd $WorkingDirectory/SNP_analysis/genomics_general
# I created this conda environment outside of the script. It includes the modules required for genomics_general
conda activate ThamnophisPopGen 
# Step 2: Create the geno files to be used by the genomics_general program popgenWindows.py
parseVCF.py -i $WorkingDirectory/variantFiltration/Full_Exons.vcf.gz | bgzip > Full_Exons.geno.gz
parseVCF.py -i $WorkingDirectory/variantFiltration/Full_CDS_missense.vcf.gz | bgzip > Full_CDS_missense.geno.gz
# Step 3: Create windows file (this is already done for the genes, just needs to be done for missense SNPs)
cd $WorkingDirectory/variantFiltration
gunzip Full_CDS_missense.vcf.gz
vcf2bed < Full_CDS_missense.vcf > $WorkingDirectory/References/Full_CDS_missense.bed
bgzip -f Full_CDS_missense.vcf
bcftools index -f Full_CDS_missense.vcf.gz
# Step 3: Create variables to be used for popGenWindows.py
cd $WorkingDirectory/SNP_analysis/genomics_general
pops=$WorkingDirectory/SNP_analysis/Populations/samples-by-population.txt
# COMPLETE python $pythonScripts/addEcotypes2.py $pops $WorkingDirectory/SNP_analysis/Populations/samples-by-ecotype.txt
ecos=$WorkingDirectory/SNP_analysis/Populations/samples-by-ecotype.txt
genes=$WorkingDirectory/SNP_analysis/genomics_general/Full_Exons.geno.gz
sites=$WorkingDirectory/SNP_analysis/genomics_general/Full_CDS_missense.geno.gz
bedfile_genes=$WorkingDirectory/References/SeqCap_CapturedGenes_abbrev.bed
bedfile_missense=$WorkingDirectory/References/Full_CDS_missense.bed
# Step 4: Get population genetics statistics for target genes from transcribed windows for each gene
popgenWindows.py --popsFile $pops --windCoords $bedfile_genes -g $genes -o genes.popgen.csv.gz -f phased -m 1 -T 8 --windType predefined --writeFailedWindows -p MAH -p MER -p PVM -p SUM -p STO -p CHR -p RON -p ROC -p ELF -p NAM -p MAR -p PAP
popgenWindows.py --popsFile $ecos --windCoords $bedfile_genes -g $genes -o genes.ecogen.csv.gz -f phased -m 1 -T 8 --windType predefined --writeFailedWindows -p L -p M
# Step 5: Get population genetics statistics for target genes for each missense site 
popgenWindows.py --popsFile $pops --windCoords $bedfile_missense -g $sites -o missense.popgen.csv.gz -f phased -m 1 -T 8 --windType predefined --writeFailedWindows -p MAH -p MER -p PVM -p SUM -p STO -p CHR -p RON -p ROC -p ELF -p NAM -p MAR -p PAP
popgenWindows.py --popsFile $ecos --windCoords $bedfile_missense -g $sites -o missense.ecogen.csv.gz -f phased -m 1 -T 8 --windType predefined --writeFailedWindows -p L -p M
# Step 6: 
echo "geneID" > Header.txt
awk '{print $4}' $WorkingDirectory/References/SeqCap_CapturedGenes.bed > GeneIDs.txt
cat Header.txt GeneIDs.txt > temp.txt ; mv temp.txt GeneIDs.txt
gunzip genes.popgen.csv.gz
gunzip genes.ecogen.csv.gz
gunzip missense.popgen.csv.gz
gunzip missense.ecogen.csv
genomics_general_add-genes GeneIDs.txt genes.popgen.csv
genomics_general_add-genes GeneIDs.txt genes.ecogen.csv
genomics_general_add-genes GeneIDs.txt missense.popgen.csv
genomics_general_add-genes GeneIDs.txt missense.ecogen.csv
}

function genomics_general_add-genes {
paste $1 $2 > temp.csv ; mv temp.csv $2
sed -i "s/\t/,/" $2
head -n -4 $2 > temp.txt ; mv temp.txt $2
}

function pixyPopGen {
  # mkdir -p $WorkingDirectory/SNP_analysis/Populations/pixyPops/$1/$2
  # cp $WorkingDirectory/SNP_analysis/Populations/pairwisePops/$1/exons/$2/"$2"_Exons_"$1".vcf* \
  #    $WorkingDirectory/SNP_analysis/Populations/pixyPops/$1/$2/
  cp $WorkingDirectory/SNP_analysis/Populations/pairwisePops/$2/missense/$1/"$1"_CDS_"$2".vcf* \
     $WorkingDirectory/SNP_analysis/Populations/pixyPops/$2/$1/
  cd $WorkingDirectory/SNP_analysis/Populations/pixyPops/$2/$1/
  mv ../../"$2".pix ../
  rm pixy_tmpfile*
  bgzip "$1"_CDS_"$2".vcf
  tabix -p vcf "$1"_CDS_"$2".vcf.gz
  bgzip "$1"_Exons_"$2".vcf
  tabix -p vcf "$1"_Exons_"$2".vcf.gz
  conda activate pixy
  pixy --stats pi fst dxy \
      --vcf "$1"_Exons_"$2".vcf.gz \
      --populations ../"$2".pix \
      --bed_file $WorkingDirectory/References/GeneBEDs/"$1"_Exons.bed \
      --output_prefix "$1"_"$2" \
      --n_core 1 \
      --bypass_invariant_check 'yes'
  conda deactivate
}

function getTranscriptLengths {
echo $1
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
if [[ "$1" == "CDS" ]]; then
  python $pythonScripts/getCDSlength.py $WorkingDirectory/References/SeqCap_Captured"$1".gff Full_"$1""$2"_Nvariants.txt Full_"$1""$2"_TranscriptLengths.txt Log.txt
else
  python $pythonScripts/getTranscriptLengths.py $WorkingDirectory/References/SeqCap_Captured"$1".gff Full_"$1""$2"_Nvariants.txt Full_"$1""$2"_TranscriptLengths.txt
fi
python $pythonScripts/getVariableRegionsGFF.py Full_"$1""$2"_TranscriptLengths.txt $WorkingDirectory/References/SeqCap_Captured"$1".gff $WorkingDirectory/References/Full_Variable"$1""$2".gff Full_"$1""$2"_variableGenes.txt

}

function vcf2faa {
# This function takes a gene-specific vcf with multiple samples and turns it into fasta files (both fna and faa) for each sample
# Prepare GFF for gffread (the gff must be in a particular format in order to work in gffread)
# DONE python $pythonScripts/modifyGFF_gffread.py $WorkingDirectory/References/Full_VariableCDS.gff $WorkingDirectory/References/Full_VariableCDS_gffread.gff
mkdir -p $WorkingDirectory/SNP_analysis/vcf2fasta
cp $WorkingDirectory/variantFiltration/Full_CDS.vcf.gz $WorkingDirectory/SNP_analysis/vcf2fasta
cd $WorkingDirectory/SNP_analysis/vcf2fasta
gunzip Full_CDS.vcf.gz
cp $WorkingDirectory/variantFiltration/Full_CDS.vcf.gz .
bcftools index -f Full_CDS.vcf.gz
bcftools query -l Full_CDS.vcf.gz > sampleList.txt
# Add the reference to the sample list
echo "RefSeq" >> sampleList.txt
# Loop through the sample list
while read sample
do
  cd $WorkingDirectory/SNP_analysis/vcf2fasta
  mkdir -p "$sample"
  # Create a .vcf for the sample
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" SelectVariants \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -V Full_CDS.vcf \
    -O "$sample"/"$sample".vcf \
    -sn "$sample"
  # Insert SNPs into the reference genome (this outputs your initial fasta file, which is the size of the genome and includes sites with low mapping depth)
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" FastaAlternateReferenceMaker \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -V "$sample"/"$sample".vcf \
    -O "$sample"/"$sample"_wholeGenome_wrongHeaders.fasta \
    --use-iupac-sample "$sample"
  # Change the fasta headers (to work in downstream programs)
  python $pythonScripts/changeGenomeHeaders.py $WorkingDirectory/References/TelagGenome.fasta "$sample"/"$sample"_wholeGenome_wrongHeaders.fasta "$sample"/"$sample"_wholeGenome.fasta
  # Use bedtools to mask regions with low mapping coverage
  bedtools genomecov \
    -ibam $WorkingDirectory/mappedReadsAll/"$sample".bam -bga | \
    awk '$4<2' | \
    bedtools maskfasta -fi "$sample"/"$sample"_wholeGenome.fasta -bed - -fo "$sample"/"$sample"_maskedGenome_wrongHeaders.fasta
  # Change the fasta headers again (they were modified by bedtools)
  python $pythonScripts/changeGenomeHeaders.py $WorkingDirectory/References/TelagGenome.fasta "$sample"/"$sample"_maskedGenome_wrongHeaders.fasta "$sample"/"$sample"_maskedGenome.fasta
  #+ this list used to be called "Log.txt"
  echo "$sample" >> vcf2faa_log.txt
  # Reduce the fasta to include only the targeted regions (-x is the outfile, -g is the infile, the last line is the reference)
  gffread \
    -x "$sample"/"$sample"_maskedCDS.fasta \
    -g "$sample"/"$sample"_maskedGenome.fasta \
    $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff 
  mkdir "$sample"/Sequences
  cd "$sample"/Sequences
  # Translate the fasta file to get the peptide sequence
  python $pythonScripts/parseAndTranslate.py ../"$sample"_maskedCDS.fasta "$sample"
done<sampleList.txt
# ^^ I created this list of samples from Full_CDS.vcf, using only the samples from Seq Cap or RNA seq
}

function reference2faa {
cd $WorkingDirectory/SNP_analysis/vcf2fasta
mkdir -p RefSeq
gffread \
  -x RefSeq/RefSeq_maskedCDS.fasta \
  -g $WorkingDirectory/References/TelagGenome.fasta \
  $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff 
mkdir RefSeq/Sequences
cd RefSeq/Sequences
python $pythonScripts/parseAndTranslate.py ../RefSeq_maskedCDS.fasta RefSeq
}

function moveCapturedGenes {
  cd $WorkingDirectory/SNP_analysis/vcf2fasta/RefSeq/Sequences
  mkdir -p CapturedGenes
  while read i
  do
    mv *"$i"* CapturedGenes/ 
  done<$WorkingDirectory/References/CapturedGenes.txt
}

function createMSA {
  cd $WorkingDirectory/SNP_analysis/vcf2fasta/RefSeq/Sequences/CapturedGenes
  ls *."$1" | cut -d "_" -f 2,3,4 | sort > "$2"List.txt
  mkdir -p $WorkingDirectory/SNP_analysis/"$4"/Captured"$1"
  cd $WorkingDirectory/SNP_analysis/"$4"/Captured"$1"
  while read fasta
  do
    while read sample
    do
      sed -i "s/>"$sample"_"$sample"_"$sample"_/"$sample"_/" $WorkingDirectory/SNP_analysis/vcf2fasta/"$sample"/"$3"/"$sample"_"$fasta"
      sed -i "s/>"$sample"_"$sample"_/"$sample"_/" $WorkingDirectory/SNP_analysis/vcf2fasta/"$sample"/"$3"/"$sample"_"$fasta"
      sed -i "s/>rna/>"$sample"_/" $WorkingDirectory/SNP_analysis/vcf2fasta/"$sample"/"$3"/"$sample"_"$fasta"
      cat $WorkingDirectory/SNP_analysis/vcf2fasta/"$sample"/"$3"/"$sample"_"$fasta" >> Alignment_"$fasta"
      echo "" >> Alignment_"$fasta"
    done<$WorkingDirectory/SNP_analysis/vcf2fasta/sampleList.txt
  done<$WorkingDirectory/SNP_analysis/vcf2fasta/RefSeq/"$3"/CapturedGenes/"$2"List.txt
}
   
function vcf2faa_unmasked {
# I already copied the RefSeq Sequences directory to UnmaskedSequences
cd $WorkingDirectory/SNP_analysis/vcf2fasta
for sample in `bcftools query -l SeqCap_CDS.vcf`
do
  cd $WorkingDirectory/SNP_analysis/vcf2fasta
  gffread \
    -x "$sample"/"$sample"_unmaskedCDS.fasta \
    -g "$sample"/"$sample"_wholeGenome.fasta \
    $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff 
  mkdir "$sample"/UnmaskedSequences
  cd "$sample"/UnmaskedSequences
  python $pythonScripts/parseAndTranslate.py ../"$sample"_unmaskedCDS.fasta "$sample"
done
}

