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

}

function sortVariants {
cd $WorkingDirectory/variantFiltration
cp SeqCap_"$1".vcf SeqCap_"$1"_original.vcf
bgzip SeqCap_"$1".vcf
bcftools index -f SeqCap_"$1".vcf.gz
bcftools norm -d snps SeqCap_"$1".vcf.gz -O v -o SeqCap_"$1"_dupsRemoved.vcf
bcftools sort SeqCap_"$1"_dupsRemoved.vcf.gz -O z -o SeqCap_"$1"_sorted.vcf.gz
mv SeqCap_"$1"_sorted.vcf.gz SeqCap_"$1".vcf.gz
rm SeqCap_"$1"_dupsRemoved.vcf.gz
bcftools index -f SeqCap_"$1".vcf.gz
# I also did this with the SeqCap_HardFilterStep3.vcf - which is the "Genes" vcf
# The file named SeqCap_Genes.vcf.gz is sorted with dups removed from the HardFilterStep3.vcf file
}

function functionalAnnotation {
#+ COMPLETED # Step 1: Download and install
#+ COMPLETED   cd ~
#+ COMPLETED   wget wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
#+ COMPLETED   gunzip snpEff_latest_core.zip
#+ COMPLETED # Step 2: Create genome annotation database
#+ COMPLETED   cd snpEff
#+ COMPLETED   mkdir data
#+ COMPLETED   cd data
#+ COMPLETED   mkdir rThaEle1 genomes
#+ COMPLETED   cp $WorkingDirectory/References/SeqCap_CapturedCDS.gff rThaEle1/genes.gff
#+ COMPLETED   cp $WorkingDirectory/References/TelagGenome.fasta genomes/rThaEle1.fa
#+ COMPLETED   cd rThaEle1
#+ COMPLETED   bgzip genes.gff
#+ COMPLETED   cd ../genomes
#+ COMPLETED   bgzip rThaEle1.gff # This could take a while, you may need to submit a job
#+ COMPLETED   cd ~/snpEff
#+ COMPLETED   java -jar /tools/snpeff-4.3p/snpEff.jar build -gff3 -v rThaEle1
#+ COMPLETED # Step 3: Run snpEff
cd $WorkingDirectory/variantFiltration
# remove singletons
#+ COMPLETED  vcftools --mac 2 --vcf Full_CDS.vcf --recode --recode-INFO-all --out Full_CDS_noSingletons.vcf
#+ COMPLETED  mv Full_CDS_noSingletons.vcf.recode.vcf Full_CDS_noSingletons.vcf
#+ COMPLETED  java -jar /tools/snpeff-4.3p/snpEff.jar -c ~/snpEff/snpEff.config -v rThaEle1 Full_CDS_noSingletons.vcf > Full_CDS_ann.vcf
#+ COMPLETED   # This didn't work with the raw gff downloaded from genbank (TelagGenome.gff). 
#+ COMPLETED   # Instead, I had to use the SeqCap_CapturedCDS.gff file I modified to only include CDS in genes of interest.
# Step 4: Pull out missense and synonymous mutations
#+ COMPLETED  awk '/^#|missense_variant/' Full_CDS_ann.vcf > Full_CDS_missense.vcf
#+ COMPLETED  awk '/^#|synonymous_variant/' Full_CDS_ann.vcf > Full_CDS_synonymous.vcf
#+ COMPLETED  echo "total CDS SNP count: $(grep -v "^#" Full_CDS_ann.vcf | wc -l)" >> Log.txt
#+ COMPLETED  echo "missense SNP count: $(grep -v "^#" Full_CDS_missense.vcf | wc -l)" >> Log.txt
#+ COMPLETED  echo "synonymous SNP count: $(grep -v "^#" Full_CDS_synonymous.vcf | wc -l)" >> Log.txt

   bgzip Full_CDS_ann.vcf
   bcftools index -f Full_CDS_ann.vcf.gz
   bcftools view -R $WorkingDirectory/References/IILS_CapturedCDS.bed.gz Full_CDS_ann.vcf.gz -O v -o Full_IILS_CDS_ann.vcf


  awk '/^#|missense_variant/' Full_IILS_CDS_ann.vcf > Full_IILS_CDS_missense.vcf
  awk '/^#|synonymous_variant/' Full_IILS_CDS_ann.vcf > Full_IILS_CDS_synonymous.vcf
  echo "IILS total CDS SNP count: $(grep -v "^#" Full_IILS_CDS_ann.vcf | wc -l)" >> Log.txt
  echo "IILS missense SNP count: $(grep -v "^#" Full_IILS_CDS_missense.vcf | wc -l)" >> Log.txt
  echo "IILS synonymous SNP count: $(grep -v "^#" Full_IILS_CDS_synonymous.vcf | wc -l)" >> Log.txt

}

function getGeneVariants {
# Get CDS SNPs and prepare for extraction
mkdir -p $WorkingDirectory/SNP_analysis/variantsByGene/"$1"
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1"
cp $WorkingDirectory/variantFiltration/SeqCap_"$1".vcf.gz .
bcftools index -f SeqCap_"$1".vcf.gz
# Create bed file for each gene
cp $WorkingDirectory/References/SeqCap_CapturedGenes.bed.gz .
gunzip SeqCap_CapturedGenes.bed.gz
python ~/SeqCap/pythonScripts/parseBED.py SeqCap_CapturedGenes.bed SeqCap_"$1"_Captures.txt "$1"
sort -u SeqCap_"$1"_Captures.txt > SeqCap_"$1"_CapturedGeneList.txt
# Extract SNPs by gene from vcf
locivar=0
while read i
  do
  mkdir -p "$i"
  mv "$i"_"$1".bed "$i"/
  bcftools view -R "$i"/"$i"_"$1".bed SeqCap_"$1".vcf.gz > "$i"/"$i"_"$1".vcf
  locusvar="$(grep -v "^#" "$i"/"$i"_"$1".vcf | wc -l)"
  echo "$i	$locusvar" >> SeqCap_"$1"_Nvariants.txt
  locivar=$((locusvar + locivar))
done<SeqCap_"$1"_CapturedGeneList.txt
  echo "total variants: $locivar" >> SeqCap_"$1"_Nvariants.txt
}

function getTranscriptLengths {
cd $WorkingDirectory/SNP_analysis/variantsByGene/"$1"
python ~/SeqCap/pythonScripts/getTranscriptLengths.py $WorkingDirectory/References/SeqCap_Captured"$1".gff SeqCap_"$1"_Nvariants.txt SeqCap_"$1"_TranscriptLengths.txt
python ~/SeqCap/pythonScripts/getVariableRegionsGFF.py SeqCap_"$1"_TranscriptLengths.txt $WorkingDirectory/References/SeqCap_Captured"$1".gff $WorkingDirectory/References/SeqCap_Variable"$1".gff SeqCap_"$1"_variableGenes.txt
}

function vcf2faa {
# Prepare GFF for gffread
python ~/SeqCap/pythonScripts/modifyGFF_gffread.py $WorkingDirectory/References/SeqCap_VariableCDS.gff $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff
mkdir -p $WorkingDirectory/SNP_analysis/vcf2fasta
cp $WorkingDirectory/variantFiltration/SeqCap_CDS.vcf.gz $WorkingDirectory/SNP_analysis/vcf2fasta
cd $WorkingDirectory/SNP_analysis/vcf2fasta
gunzip SeqCap_CDS.vcf.gz
cp $WorkingDirectory/variantFiltration/SeqCap_CDS.vcf.gz .
bcftools index -f SeqCap_CDS.vcf.gz
for sample in `bcftools query -l SeqCap_CDS.vcf`
do
  cd $WorkingDirectory/SNP_analysis/vcf2fasta
  mkdir -p "$sample"
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" SelectVariants \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -V SeqCap_CDS.vcf \
    -O "$sample"/"$sample".vcf \
    -sn "$sample"
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" FastaAlternateReferenceMaker \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -V "$sample"/"$sample".vcf \
    -O "$sample"/"$sample"_wholeGenome_wrongHeaders.fasta \
    --use-iupac-sample "$sample"
  python ~/SeqCap/pythonScripts/changeGenomeHeaders.py $WorkingDirectory/References/TelagGenome.fasta "$sample"/"$sample"_wholeGenome_wrongHeaders.fasta "$sample"/"$sample"_wholeGenome.fasta
  bedtools genomecov \
    -ibam $WorkingDirectory/mappedReadsAll/"$sample".bam -bga | \
    awk '$4<2' | \
    bedtools maskfasta -fi "$sample"/"$sample"_wholeGenome.fasta -bed - -fo "$sample"/"$sample"_maskedGenome_wrongHeaders.fasta
  python ~/SeqCap/pythonScripts/changeGenomeHeaders.py $WorkingDirectory/References/TelagGenome.fasta "$sample"/"$sample"_maskedGenome_wrongHeaders.fasta "$sample"/"$sample"_maskedGenome.fasta
  echo "$sample" >> Log.txt
  gffread \
    -x "$sample"/"$sample"_maskedCDS.fasta \
    -g "$sample"/"$sample"_maskedGenome.fasta \
    $WorkingDirectory/References/SeqCap_VariableCDS_gffread.gff 
  mkdir "$sample"/Sequences
  cd "$sample"/Sequences
  python ~/SeqCap/pythonScripts/parseAndTranslate.py ../"$sample"_maskedCDS.fasta "$sample"
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
python ~/SeqCap/pythonScripts/parseAndTranslate.py ../RefSeq_maskedCDS.fasta RefSeq
# echo RefSeq >> Log.txt
}

function createMSA {
# I changed the directory from Sequences to UnmaskedSequences, make sure to change back if you need masking
  cd $WorkingDirectory/SNP_analysis/vcf2fasta/RefSeq/"$3"
  mkdir -p CapturedGenes
  while read i
  do
    mv *"$i"* CapturedGenes/ 
  done<$WorkingDirectory/References/CapturedGenes.txt
  cd CapturedGenes
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
    done<$WorkingDirectory/SNP_analysis/vcf2fasta/Log.txt
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
  python ~/SeqCap/pythonScripts/parseAndTranslate.py ../"$sample"_unmaskedCDS.fasta "$sample"
done
}

# MAIN
loadModules
WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
#+ COMPLETED combineDatasets
#+ COMPLETED sortVariants CDS
#+ COMPLETED sortVariants Exons
#+ COMPLETED sortVariants Genes
#+ COMPLETED 
#+ COMPLETED getGeneVariants CDS
#+ COMPLETED getGeneVariants Exons
#+ COMPLETED getGeneVariants Genes
#+ COMPLETED 
#+ COMPLETED getTranscriptLengths CDS
#+ COMPLETED getTranscriptLengths Exons
#+ COMPLETED getTranscriptLengths Genes
#+ COMPLETED 
#+ COMPLETED vcf2faa
#+ COMPLETED reference2faa
# createMSA fna gene
#+ COMPLETED createMSA faa protein Sequences maskedMSA
#+ COMPLETED vcf2faa_unmasked
#+ COMPLETED createMSA faa protein UnmaskedSequences unmaskedMSA
#+ COMPLETED combineDatasets
functionalAnnotation
