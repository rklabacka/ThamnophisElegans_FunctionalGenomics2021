#!/bin/sh

# This is the bash script used for SNP analysis of the results
# from the reads2vcf.sh script

# Prepared by Randy Klabacka

# -- Job details for Hopper cluster at Auburn University -- ##
#Give job a name
#PBS -N FullScript_May2020

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=1:ppn=10,walltime=05:00:00:00
#-- Indicate if\when you want to receive email about your job
#-- The directive below sends email if the job is (a) aborted, 
#-- when it (b) begins, and when it (e) ends
#PBS -m abe rlk0015@auburn.edu

echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------

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
bgzip WGS_Genes.vcf
bgzip WGS_Exons.vcf
bgzip WGS_CDS.vcf
bcftools index -f WGS_Genes.vcf.gz
bcftools index -f WGS_Exons.vcf.gz
bcftools index -f WGS_CDS.vcf.gz
echo "merging WGS data to RNA-Seq+Seq-Cap" >> Log.txt
bcftools merge SeqCap_CDS.vcf.gz WGS_CDS.vcf.gz -O v -o Full_CDS.vcf
bcftools merge SeqCap_Exons.vcf.gz WGS_Exons.vcf.gz -O v -o Full_Exons.vcf
bcftools merge SeqCap_Genes.vcf.gz WGS_Genes.vcf.gz -O v -o Full_Genes.vcf
bgzip Full_CDS.vcf
bgzip Full_Exons.vcf
bgzip Full_Genes.vcf
bcftools index -f Full_CDS.vcf
bcftools index -f Full_Exons.vcf
bcftools index -f Full_Genes.vcf
}

function sortVariants {
cd $WorkingDirectory/variantFiltration
cp "$1".vcf.gz "$1"_original.vcf.gz
bcftools index -f "$1".vcf.gz
bcftools norm -d snps "$1".vcf.gz -O v -o "$1"_dupsRemoved.vcf
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
  mv *.txt allPops
  cd allPops
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
vcf2bed < "$1".vcf > "$1"_variants.bed
bgzip "$1"_variants.bed
bgzip "$1".vcf
}

function functionalAnnotation {
# Step 1: Download and install
  cd ~
  wget wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
  gunzip snpEff_latest_core.zip
# Step 2: Create genome annotation database
  cd snpEff
  mkdir -p data
  cd data
  mkdir -p rThaEle1 genomes
  cp $WorkingDirectory/References/"$1"_CapturedCDS.gff rThaEle1/genes.gff
  cp $WorkingDirectory/References/TelagGenome.fasta genomes/rThaEle1.fa
  cd ~/snpEff
  java -jar /tools/snpeff-4.3p/snpEff.jar build -gff3 -v rThaEle1
# Step 3: Run snpEff
  cd $WorkingDirectory/variantFiltration
  gunzip "$1"_CDS.vcf.gz
  java -jar /tools/snpeff-4.3p/snpEff.jar -c ~/snpEff/snpEff.config -v rThaEle1 "$1"_CDS.vcf > "$1"_CDS_ann.vcf
  # This didn't work with the raw gff downloaded from genbank (TelagGenome.gff). 
  # Instead, I had to use the "$1"_CapturedCDS.gff file I modified to only include CDS in genes of interest.
# Step 4: Pull out missense and synonymous mutations
  awk '/^#|missense_variant/' "$1"_CDS_ann.vcf > "$1"_CDS_missense.vcf
  awk '/^#|synonymous_variant/' "$1"_CDS_ann.vcf > "$1"_CDS_synonymous.vcf
  echo "total CDS SNP count: $(grep -v "^#" "$1"_CDS_ann.vcf | wc -l)" >> Log.txt
  echo "missense SNP count: $(grep -v "^#" "$1"_CDS_missense.vcf | wc -l)" >> Log.txt
  echo "synonymous SNP count: $(grep -v "^#" "$1"_CDS_synonymous.vcf | wc -l)" >> Log.txt
  bgzip "$1"_CDS.vcf.gz
  bgzip "$1"_CDS_missense.vcf
  bgzip "$1"_CDS_synonymous.vcf
  bcftools index -f "$1"_CDS_missense.vcf.gz
  bcftools index -f "$1"_CDS_synonymous.vcf.gz
  bgzip "$1"_CDS_ann.vcf
  bcftools index -f "$1"_CDS_ann.vcf.gz

}

function getGeneVariants {
# Get CDS SNPs and prepare for extraction
mkdir -p $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
cp $WorkingDirectory/variantFiltration/Full_"$1""$2".vcf.gz .
bcftools index -f Full_"$1""$2".vcf.gz
#+ COMPLETED # Create bed file for each gene
#+ COMPLETED getBEDbyGene $1 $2
# Extract SNPs by gene from vcf
echo "begin SNP extraction by gene from vcf" >> Log.txt
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
# WARNING: The following command has not been verified
# within the getGeneVariants function.
## -- Previously the code in the function getVCFbyGene
## -- was included as hard code within the getGeneVariants
## -- function. I created a separate function for
## -- getVCFbyGene when I needed to use the same process
## -- for the createPairwiseVCFs function. That being said,
## -- it should work fine.
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
# echo "    Entered getGeneVariants for $3"
# Extract SNPs by gene from vcf
locivar=0
while read i
  do
 #  echo "        Gene: $i"
  mkdir -p "$i"
  bcftools view -R $WorkingDirectory/References/GeneBEDs/"$i"_"$1".bed "$2" > "$i"/"$i"_"$1""$4"_"$3".vcf
  locusvar="$(grep -v "^#" "$i"/"$i"_"$1""$4"_"$3".vcf| wc -l)"
  echo "$i	$locusvar" >> Full_"$1""$4"_Nvariants.txt
  locivar=$((locusvar + locivar))
done<$WorkingDirectory/References/GeneBEDs/Full_CDS_CapturedGeneList.txt
  echo "total variants: $locivar" >> Full_"$1""$2"_Nvariants.txt
}

function getPairwisePopGen {
cd $WorkingDirectory
rm popGenLog.txt
mkdir -p $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats
cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops
echo "begin pairwise pop gen" >> Log.txt
while read i
do
  cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats
  rm "$i".txt
  echo -e "Population\tN\tMissenseSNPs\tSynonymousSNPs\tTranscriptSNPs\tTajD\tfst\tdxy\tpop1pi\tpop2pi" >> "$i".txt
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
    # pixy function to calculate average fst and dxy for exonic regions of each gene
    pixyPopGen $i $j
    cd $WorkingDirectory/SNP_analysis/Populations/pixyPops/$j/$i/
    fst="$(awk 'NR == 2 {print $6}' "$i"_"$j"_perGene_fst.txt)"
    fst=${fst:="NA"}
    echo -e "\tfst: ${fst}" >> $WorkingDirectory/popGenLog.txt
    dxy="$(awk 'NR == 2 {print $6}' "$i"_"$j"_perGene_dxy.txt)"
    dxy=${dxy:="NA"}
    echo -e "\tdxy: ${dxy}" >> $WorkingDirectory/popGenLog.txt
    # the pi values calculated from pixy are average pairwise differences within each population
    pop1pi="$(awk 'NR == 2 {print $5}' "$i"_"$j"_perGene_pi.txt)"
    pop1pi=${pop1pi:="NA"}
    echo -e "\tpop1pi: ${pop1pi}" >> $WorkingDirectory/popGenLog.txt
    pop2pi="$(awk 'NR == 3 {print $5}' "$i"_"$j"_perGene_pi.txt)"
    pop2pi=${pop2pi:="NA"}
    echo -e "\tpop2pi: ${pop2pi}" >> $WorkingDirectory/popGenLog.txt
    cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/"$j"/exons/"$i"
    echo -e "$j\t$n\t$misSNPcount\t$synSNPcount\t$transcriptSNPcount\t$tajD\t$fst\t$dxy\t$pop1pi\t$pop2pi" >> $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats/"$i".txt
  done<$WorkingDirectory/SNP_analysis/Populations/pairwisePops/PairwisePopsList
  cd $WorkingDirectory/SNP_analysis/Populations/pairwisePops/PopGenStats
  python $pythonScripts/addEcotypes.py "$i".txt "$i"_withEcotypes.txt
done<$WorkingDirectory/References/GeneBEDs/Full_CDS_CapturedGeneList.txt
}

function pairwisePopGen-2 {
cd $WorkingDirectory/variantFiltration
python3 ~/myTools/genomics_general/VCF_processing/parseVCF.py -i Full_Exons.vcf.gz | bgzip > Full_Exons.geno.gz


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
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
python $pythonScripts/getTranscriptLengths.py $WorkingDirectory/References/SeqCap_Captured"$1".gff Full_"$1""$2"_Nvariants.txt Full_"$1""$2"_TranscriptLengths.txt
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

