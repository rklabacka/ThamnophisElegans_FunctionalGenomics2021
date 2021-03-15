#!/bin/sh

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

function loadModules {
  module load gnu-parallel/20160322 
  module load samtools/1.3.1
  module load xz/5.2.2
  module load htslib/1.3.3
  module load python/3.6.4
  module load gatk/4.1.7.0
  module load R/3.6.0
  module load bedtools/2.29.0
  module load bcftools/1.3.2
  module load vcftools/v0.1.17
  module load gffread/2
}

function combineDatasets {

cd $WorkingDirectory/variantFiltration
# Jessica's WGS variants at variable sites discovered from q.FullScript_May2020.sh WGS_Genes.recode.vcf.gz WGS_Exons.recode.vcf.gz and WGS_CDS.recode.vcf.gz were copied here frome box already
bcftools index -f SeqCap_CDS.vcf.gz
bcftools index -f SeqCap_Exons.vcf.gz
bcftools index -f SeqCap_Genes.vcf.gz
bcftools index -f CDS_WGS.recode.vcf.gz
bcftools index -f Exons_WGS.recode.vcf.gz
bcftools index -f Genes_WGS.recode.vcf.gz
bcftools merge SeqCap_CDS.vcf.gz CDS_WGS.recode.vcf.gz -O v -o Full_CDS.vcf
bcftools merge SeqCap_Exons.vcf.gz Exons_WGS.recode.vcf.gz -O v -o Full_Exons.vcf
bcftools merge SeqCap_Genes.vcf.gz Genes_WGS.recode.vcf.gz -O v -o Full_Genes.vcf
bgzip Full_CDS.vcf
bgzip Full_Exons.vcf
bgzip Full_Genes.vcf
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
  for sample in `bcftools query -l $WorkingDirectory/variantFiltration/Full_CDS.vcf.gz`
  do
    echo "$sample" >> Samples
  done
  python $pythonScripts/parsePopulations.py Samples.txt
  rm *PIK*.txt
  mkdir -p pairwisePops
  mkdir -p allPops
  mv *_*.txt pairwisePops
  mv *.txt allPops
} 
function getVariantBED {
cd $WorkingDirectory/variantFiltration
gunzip "$1".vcf.gz
vcf2bed < "$1".vcf > "$1"_variants.bed
bgzip "$1"_variants.bed
bgzip "$1".vcf
}

function functionalAnnotation {
#+ COMPLETED # Step 1: Download and install
  cd ~
#+ COMPLETED   wget wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
#+ COMPLETED   gunzip snpEff_latest_core.zip
#+ COMPLETED # Step 2: Create genome annotation database
#+ COMPLETED   cd snpEff
#+ COMPLETED   mkdir -p data
#+ COMPLETED   cd data
#+ COMPLETED   mkdir -p rThaEle1 genomes
#+ COMPLETED   cp $WorkingDirectory/References/"$1"_CapturedCDS.gff rThaEle1/genes.gff
#+ COMPLETED #+ COMPLETED cp $WorkingDirectory/References/TelagGenome.fasta genomes/rThaEle1.fa
#+ COMPLETED   cd ~/snpEff
#+ COMPLETED   java -jar /tools/snpeff-4.3p/snpEff.jar build -gff3 -v rThaEle1
# Step 3: Run snpEff
cd $WorkingDirectory/variantFiltration
# remove singletons
#+ COMPLETED  vcftools --mac 2 --vcf "$1"_CDS.vcf --recode --recode-INFO-all --out "$1"_CDS_noSingletons.vcf
#+ COMPLETED  mv "$1"_CDS_noSingletons.vcf.recode.vcf "$1"_CDS_noSingletons.vcf
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
#   bgzip "$1"_CDS_ann.vcf
#   bcftools index -f "$1"_CDS_ann.vcf.gz
#   bcftools view -R $WorkingDirectory/References/IILS_CapturedCDS.bed.gz "$1"_CDS_ann.vcf.gz -O v -o "$1"_IILS_CDS_ann.vcf
#
#
#  awk '/^#|missense_variant/' "$1"_IILS_CDS_ann.vcf > "$1"_IILS_CDS_missense.vcf
#  awk '/^#|synonymous_variant/' "$1"_IILS_CDS_ann.vcf > "$1"_IILS_CDS_synonymous.vcf
#  echo "IILS total CDS SNP count: $(grep -v "^#" "$1"_IILS_CDS_ann.vcf | wc -l)" >> Log.txt
#  echo "IILS missense SNP count: $(grep -v "^#" "$1"_IILS_CDS_missense.vcf | wc -l)" >> Log.txt
#  echo "IILS synonymous SNP count: $(grep -v "^#" "$1"_IILS_CDS_synonymous.vcf | wc -l)" >> Log.txt

}

function getGeneVariants {
# Get CDS SNPs and prepare for extraction
mkdir -p $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
cp $WorkingDirectory/variantFiltration/Full_"$1""$2".vcf.gz .
bcftools index -f Full_"$1""$2".vcf.gz
# Create bed file for each gene
cp $WorkingDirectory/References/SeqCap_CapturedGenes.bed.gz .
gunzip SeqCap_CapturedGenes.bed.gz
python $pythonScripts/parseBED.py SeqCap_CapturedGenes.bed Full_"$1""$2"_Captures.txt "$1"
sort -u Full_"$1""$2"_Captures.txt > Full_"$1""$2"_CapturedGeneList.txt
# Extract SNPs by gene from vcf
locivar=0
while read i
  do
  mkdir -p "$i"
  mv "$i"_"$1".bed "$i"/
  bcftools view -R "$i"/"$i"_"$1".bed Full_"$1""$2".vcf.gz > "$i"/"$i"_"$1""$2".vcf
  locusvar="$(grep -v "^#" "$i"/"$i"_"$1""$2".vcf | wc -l)"
  echo "$i	$locusvar" >> Full_"$1""$2"_Nvariants.txt
  locivar=$((locusvar + locivar))
done<Full_"$1""$2"_CapturedGeneList.txt
  echo "total variants: $locivar" >> Full_"$1"_Nvariants.txt
}

function getTranscriptLengths {
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1""$2"
python $pythonScripts/getTranscriptLengths.py $WorkingDirectory/References/SeqCap_Captured"$1".gff Full_"$1""$2"_Nvariants.txt Full_"$1""$2"_TranscriptLengths.txt
python $pythonScripts/getVariableRegionsGFF.py Full_"$1""$2"_TranscriptLengths.txt $WorkingDirectory/References/SeqCap_Captured"$1".gff $WorkingDirectory/References/Full_Variable"$1""$2".gff Full_"$1""$2"_variableGenes.txt
}

function vcf2faa {
# Prepare GFF for gffread
python $pythonScripts/modifyGFF_gffread.py $WorkingDirectory/References/Full_VariableCDS.gff $WorkingDirectory/References/Full_VariableCDS_gffread.gff
mkdir -p $WorkingDirectory/SNP_analysis/vcf2fasta
cp $WorkingDirectory/variantFiltration/Full_CDS.vcf.gz $WorkingDirectory/SNP_analysis/vcf2fasta
cp $WorkingDirectory/variantFiltration/SeqCap_CDS.vcf.gz $WorkingDirectory/SNP_analysis/vcf2fasta
cd $WorkingDirectory/SNP_analysis/vcf2fasta
gunzip Full_CDS.vcf.gz
gunzip SeqCap_CDS.vcf.gz
cp $WorkingDirectory/variantFiltration/Full_CDS.vcf.gz .
cp $WorkingDirectory/variantFiltration/SeqCap_CDS.vcf.gz .
bcftools index -f Full_CDS.vcf.gz
bcftools index -f SeqCap_CDS.vcf.gz
#+ COMPLETED echo "RefSeq" >> SampleList.txt
for sample in `bcftools query -l SeqCap_CDS.vcf`
do
  cd $WorkingDirectory/SNP_analysis/vcf2fasta
  mkdir -p "$sample"
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" SelectVariants \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -V Full_CDS.vcf \
    -O "$sample"/"$sample".vcf \
    -sn "$sample"
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" FastaAlternateReferenceMaker \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -V "$sample"/"$sample".vcf \
    -O "$sample"/"$sample"_wholeGenome_wrongHeaders.fasta \
    --use-iupac-sample "$sample"
  python $pythonScripts/changeGenomeHeaders.py $WorkingDirectory/References/TelagGenome.fasta "$sample"/"$sample"_wholeGenome_wrongHeaders.fasta "$sample"/"$sample"_wholeGenome.fasta
  bedtools genomecov \
    -ibam $WorkingDirectory/mappedReadsAll/"$sample".bam -bga | \
    awk '$4<2' | \
    bedtools maskfasta -fi "$sample"/"$sample"_wholeGenome.fasta -bed - -fo "$sample"/"$sample"_maskedGenome_wrongHeaders.fasta
  python $pythonScripts/changeGenomeHeaders.py $WorkingDirectory/References/TelagGenome.fasta "$sample"/"$sample"_maskedGenome_wrongHeaders.fasta "$sample"/"$sample"_maskedGenome.fasta
  #+ COMPLETED This is complete- also I noticed that the SeqCap_CDS.vcf file I used for this command contains duplicate samples that I needed to delete manually
  #+ COMPLETED Also, this list used to be called "Log.txt"
  #+ COMPLETED echo "$sample" >> SampleList.txt
  gffread \
    -x "$sample"/"$sample"_maskedCDS.fasta \
    -g "$sample"/"$sample"_maskedGenome.fasta \
    $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff 
  mkdir "$sample"/Sequences
  cd "$sample"/Sequences
  python $pythonScripts/parseAndTranslate.py ../"$sample"_maskedCDS.fasta "$sample"
done
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

function reference2faa {
cd $WorkingDirectory/SNP_analysis/vcf2fasta
echo "Outgroup" >> SampleList.txt
mkdir -p "$1"
gffread \
  -x "$1"/"$1"_maskedCDS.fasta \
  -g $WorkingDirectory/References/TelagGenome.fasta \
  $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff 
mkdir "$1"/Sequences
cd "$1"/Sequences
python $pythonScripts/parseAndTranslate.py ../"$1"_maskedCDS.fasta "$1"
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
    done<$WorkingDirectory/SNP_analysis/vcf2fasta/SampleList.txt
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

# MAIN
loadModules
WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
#+ COMPLETED combineDatasets
#+ COMPLETED sortVariants SeqCap_CDS   SeqCap_IndividualsToRemove 
#+ COMPLETED sortVariants SeqCap_Exons  SeqCap_IndividualsToRemoves
#+ COMPLETED sortVariants SeqCap_Genes  SeqCap_IndividualsToRemoves
#+ COMPLETED sortVariants Full_CDS Full_IndividualsToRemove
#+ COMPLETED sortVariants Full_Exons Full_IndividualsToRemove
#+ COMPLETED sortVariants Full_Genes Full_IndividualsToRemove

#+ COMPLETED getVariantBED SeqCap_CDS
#+ COMPLETED getVariantBED SeqCap_Exons
#+ COMPLETED getVariantBED SeqCap_Genes
#+ 
#+ COMPLETED and deleted... functionalAnnotation SeqCap
#+ COMPLETED functionalAnnotation Full
#+
#+ COMPLETED getGeneVariants CDS
#+ COMPLETED getGeneVariants Exons
#+ COMPLETED getGeneVariants Genes
#+ COMPLETED getGeneVariants CDS _missense
#+ COMPLETED getGeneVariants CDS _synonymous
#+ 
#+ COMPLETED getTranscriptLengths CDS
#+ COMPLETED getTranscriptLengths Exons
#+ COMPLETED getTranscriptLengths Genes
#+ COMPLETED getTranscriptLengths CDS _missense
#+ COMPLETED getTranscriptLengths CDS _synonymous
#+ 
#+ COMPLETED vcf2faa
#+ COMPLETED reference2faa
#+ COMPLETED moveCapturedGenes
#+ COMPLETED createMSA faa protein Sequences maskedMSA
createMSA fna transcript Sequences maskedMSA
#+ COMPLETED vcf2faa_unmasked
#+ COMPLETED createMSA faa protein UnmaskedSequences unmaskedMSA
#+ COMPLETED combineDatasets
