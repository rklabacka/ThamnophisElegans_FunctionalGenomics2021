#!/bin/sh

# This is the bash script used to process raw reads for the garter snake 
# targeted gene network study (for gene expression and populationi genetics)

# Prepared by Randy Klabacka

# -- Job Details for Hopper cluster at Auburn University-- #
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

# -- Functions -- #
## --- To see implementation of functions see q.main.sh--- ###

function createWorkingEnvironment {
  #create working directory ###
    mkdir -p /scratch/rlk0015/Telag/May2020/WorkingDirectory 
  # DNA
    mkdir -p $WorkingDirectory/rawReadsDNA
    mkdir -p $WorkingDirectory/cleanReadsDNA
    mkdir -p $WorkingDirectory/mappedReadsDNA
    mkdir -p $WorkingDirectory/StatsDNA
    mkdir -p $WorkingDirectory/GATKDNA
    mkdir -p $WorkingDirectory/SNPTablesDNA
    mkdir -p $WorkingDirectory/SequenceTablesDNA
  # RNA
    mkdir -p $WorkingDirectory/rawReadsRNA
    mkdir -p $WorkingDirectory/cleanReadsRNA
    mkdir -p $WorkingDirectory/mappedReadsRNA
    mkdir -p $WorkingDirectory/StatsRNA
    mkdir -p $WorkingDirectory/GATKRNA
    mkdir -p $WorkingDirectory/SNPTablesRNA
    mkdir -p $WorkingDirectory/SequenceTablesRNA
  # Joint
    mkdir -p $WorkingDirectory/mappedReadsAll
  # Utilities
    mkdir -p $WorkingDirectory/References
    echo "RLK_report: directory created: $WorkingDirectory with rawReads and cleanReads sub directories"
}
   
function copyRawReadsDNA {
  ### copy over raw reads ### 
    cd /home/shared/tss0019_lab/SeqCap_GarterSnake2012/
    for i in {1..96}
    do
      cp Sample_HTAdapter"$i"/*.fastq.gz $WorkingDirectory/rawReadsDNA
    done
    echo "RLK_report: SEQ-CAP RAW READ COPY COMPLETE"
  ### copy over adapter file ###
   cp /home/rlk0015/SeqCap/code/References/adapters.fa $WorkingDirectory/rawReadsDNA
   echo "RLK_report: ADAPTERS COPY COMPLETE"
  
}

function performFASTQC {
  cd $WorkingDirectory/"$1"
  ls *.fastq.gz | time parallel -j+0 --eta 'fastqc {}'
  multiqc .
  echo "RLK_report: QUALITY CHECK COMPLETE"
}
   
function performTrimmingPE {   
   cd $WorkingDirectory/rawReads"$1"
   ### paired-end trimming ###
   ls | grep "fastq.gz" | cut -d "_" -f 1,2 | sort | uniq > PE_TrimmList
   while read i
   do
   java -jar /tools/trimmomatic-0.36/trimmomatic-0.36.jar \
     PE \
     -threads 6 \
     -phred33 \
     $WorkingDirectory/rawReads"$1"/"$i"*_R1*.fastq.gz \
     $WorkingDirectory/rawReads"$1"/"$i"*_R2*.fastq.gz \
     $WorkingDirectory/cleanReads"$1"/"$i"_R1_paired.fastq.gz \ 
     $WorkingDirectory/cleanReads"$1"/"$i"_R1_unpaired.fastq.gz \
     $WorkingDirectory/cleanReads"$1"/"$i"_R2_paired.fastq.gz \
     $WorkingDirectory/cleanReads"$1"/"$i"_R2_unpaired.fastq.gz \
     ILLUMINACLIP:adapters.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:6:30 MINLEN:36
   done<PE_TrimmList
}
   
function performTrimmingSE {
  cd $WorkingDirectory/rawReads"$1"
    java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar SE -phred33 "$i"_1.fastq.gz  /scratch/GarterSnake/RNAseq_2008Data/CleanedData/"$i"_cleaned.fastq.gz  LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:36
}

function copyRef {
   cp /home/rlk0015/SeqCap/code/References/T_elegans_genome/latest_assembly_versions/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.fna.gz $WorkingDirectory/References/TelagGenome.fasta.gz
   cp /home/rlk0015/SeqCap/code/References/T_elegans_genome/latest_assembly_versions/GCF_009769535.1_rThaEle1.pri/GCF_009769535.1_rThaEle1.pri_genomic.gff.gz $WorkingDirectory/References/TelagGenome.gff.gz
   cd $WorkingDirectory/References
   gunzip TelagGenome.fasta.gz
   index T. elegans genome
   cd $WorkingDirectory/References
   bwa index -p TelagGenome -a bwtsw TelagGenome.fasta
   gunzip TelagGenome.gff.gz
}

function mapReadsDNA {
   # create list with each paired individual
   cd $WorkingDirectory/cleanReadsDNA
   ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > $WorkingDirectory/mappedReadsDNA/pairedMapList
   echo "pairedMapList created"
   cd $WorkingDirectory/mappedReadsDNA
   # while loop through the names in pairedMapList
   while read i
     do
     ### map to T. elegans genome ###
     bwa mem -t 4 -M $WorkingDirectory/References/TelagGenome \
     	$WorkingDirectory/cleanReadsDNA/"$i"_*R1_paired.fastq.gz \
     	$WorkingDirectory/cleanReadsDNA/"$i"_*R2_paired.fastq.gz > \
     	$WorkingDirectory/mappedReadsDNA/"$i"_mapped.sam
     echo "mapped DNA reads to genome"
     ###  convert .sam to .bam & sort  ###
     samtools view -@ 2 -bS $WorkingDirectory/mappedReads"$1"/"$i"*_mapped.sam | samtools sort -@ 2 -o $WorkingDirectory/mappedReads"$1"/"$i"_sorted.bam
   done<$WorkingDirectory/mappedReadsDNA/pairedMapList
}

function mapSEReadsRNA {
  # This section was performed by TSS in a separate script, and is included here only for parameter reference
  ###  Prepare the reference Index for HiSat2
    cd $WorkingDirectory/References
    gffread $WorkingDirectory/References/TelegGenome.gff -T -o $WorkingDirectory/References/TelagGenome.gtf
    extract_splice_sites.py $WorkingDirectory/References/TelagGenome.gtf > $WorkingDirectory/References/TelagGenome.ss
    extract_exons.py $WorkingDirectory/References/TelagGenome.gtf > $WorkingDirectory/References/TelagGenome.exon
  ###Create a HISAT2 index
    hisat2-build --ss $WorkingDirectory/References/.ss --exon $WorkingDirectory/References/TelagGenome.exon $WorkingDirectory/References/TelagGenome.fasta TelagGenome_index
  # Move to the 2008 SE data directory
    $WorkingDirectory/cleanReadsRNA
    wc -l *.fastq  >>  LineCount_Cleaned.txt
  ## Create list of fastq files to map
    ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list    #should list SU_1703_USR16088798L_HHHKJBBXX_L7
    cd $WorkingDirectory/mappedReadsRNA
  # copy the list of unique ids from the original files to map
  # This is from the 2008 SE directory
    cp $WorkingDirectory/cleanReadsRNA/list . 
    mkdir ballgown
  while read i;
    do
    ##HiSat  -p indicates 20 processors, --dta reports alignments for StringTieq
    hisat2 -p 20 --dta -x $WorkingDirectory/References/TelagGenome_index -U $WorkingDirectory/cleanReadsRNA/"$i"_cleaned.fastq -S "$i".sam
      ###view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
      ###sort: convert the BAM file to a sorted BAM file.
      # Example Input: HS06_GATCAG_All.sam; Output: HS06_GATCAG_sorted.bam
    samtools view -@ 20 -bS "$i".sam  | samtools sort -@ 20 -o  "$i"_sorted.bam   
    mkdir ballgown/"$i"
    # Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
    #stringtie -p 20 --rf -G "$REFDIR"/Daphnia_pulex.gtf -o "$i".gtf  -l "$i"  "$i"_sorted.bam 
        # eAB: This will run stringtie once and hopefully ONLY use the Ref annotation
    stringtie -p 20 -e -B -G $WorkingDirectory/References/TelagGenome.gtf -o ballgown/"$i"/"$i".gtf  -l "$i"  $WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam
  done<list
}

function mapPEReadsRNA {
  # This section was performed by TSS in a separate script, and is included here only for parameter reference
  while read i;
    do
    	##HiSat  -p indicates 20 processors, --dta reports alignments for StringTie --rf is the read orientation (from 
    	#mike Crowley our reads are rf since made with Agilent Stranded RNA prep)
    hisat2 -p 20 --dta --phred33 -x $WorkingDirectory/References/TelagGenome_index -1 $WorkingDirectory/cleanReadsRNA/"$i"_1_paired.fastq -2 $WorkingDirectory/cleanReadsRNA/"$i"_2_paired.fastq -S "$i".sam
        ###view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
        ###sort: convert the BAM file to a sorted BAM file.
        # Example Input: HS06_GATCAG_All.sam; Output: HS06_GATCAG_sorted.bam
    samtools view -@ 20 -bS "$i".sam  | samtools sort -@ 20 -o  "$i"_sorted.bam   
    mkdir ballgown/"$i"
      # Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
      #stringtie -p 20 --rf -G "$REFDIR"/Daphnia_pulex.gtf -o "$i".gtf  -l "$i"  "$i"_sorted.bam 
    	# eAB: This will run stringtie once and hopefully ONLY use the Ref annotation
    stringtie -p 20 -e -B -G $WorkingDirectory/References/TelagGenome.gtf -o ballgown/"$i"/"$i".gtf  -l "$i"  $WorkingDirectory/mappedReadsRNA/"$i"_sorted.bam
  done<list
}
function changeSeqCapNames {
  cd $WorkingDirectory/mappedReadsDNA
  mv HTAdapter1_sorted.bam ELF01_sorted.bam
  mv HTAdapter2_sorted.bam ELF02_sorted.bam
  mv HTAdapter3_sorted.bam ELF03_sorted.bam
  mv HTAdapter4_sorted.bam ELF04_sorted.bam
  mv HTAdapter5_sorted.bam ELF05_sorted.bam
  mv HTAdapter6_sorted.bam ELF06_sorted.bam
  mv HTAdapter7_sorted.bam ELF07_sorted.bam
  mv HTAdapter8_sorted.bam ELF08_sorted.bam
  mv HTAdapter9_sorted.bam ELF09_sorted.bam
  mv HTAdapter10_sorted.bam ELF10_sorted.bam
  mv HTAdapter11_sorted.bam ELF11_sorted.bam
  mv HTAdapter12_sorted.bam ELF12_sorted.bam
  mv HTAdapter13_sorted.bam ELF13_sorted.bam
  mv HTAdapter14_sorted.bam ELF14_sorted.bam
  mv HTAdapter15_sorted.bam ELF15_sorted.bam
  mv HTAdapter16_sorted.bam ELF16_sorted.bam
  mv HTAdapter17_sorted.bam MAH01_sorted.bam
  mv HTAdapter18_sorted.bam MAH02_sorted.bam
  mv HTAdapter19_sorted.bam MAH03_sorted.bam
  mv HTAdapter20_sorted.bam MAH04_sorted.bam
  mv HTAdapter21_sorted.bam MAH05_sorted.bam
  mv HTAdapter22_sorted.bam MAH06_sorted.bam
  mv HTAdapter23_sorted.bam MAH07_sorted.bam
  mv HTAdapter24_sorted.bam MAH08_sorted.bam
  mv HTAdapter25_sorted.bam MAH09_sorted.bam
  mv HTAdapter26_sorted.bam MAH10_sorted.bam
  mv HTAdapter27_sorted.bam MAH11_sorted.bam
  mv HTAdapter28_sorted.bam MAH12_sorted.bam
  mv HTAdapter29_sorted.bam MAH13_sorted.bam
  mv HTAdapter30_sorted.bam MAH14_sorted.bam
  mv HTAdapter31_sorted.bam MAH15_sorted.bam
  mv HTAdapter32_sorted.bam MAR01_sorted.bam
  mv HTAdapter33_sorted.bam MAR02_sorted.bam
  mv HTAdapter34_sorted.bam MAR03_sorted.bam
  mv HTAdapter35_sorted.bam MAR04_sorted.bam
  mv HTAdapter36_sorted.bam MAR05_sorted.bam
  mv HTAdapter37_sorted.bam MAR06_sorted.bam
  mv HTAdapter38_sorted.bam MAR07_sorted.bam
  mv HTAdapter39_sorted.bam MAR08_sorted.bam
  mv HTAdapter40_sorted.bam MAR09_sorted.bam
  mv HTAdapter41_sorted.bam MAR10_sorted.bam
  mv HTAdapter42_sorted.bam PAP01_sorted.bam
  mv HTAdapter43_sorted.bam PAP02_sorted.bam
  mv HTAdapter44_sorted.bam PAP03_sorted.bam
  mv HTAdapter45_sorted.bam PAP04_sorted.bam
  mv HTAdapter46_sorted.bam PAP05_sorted.bam
  mv HTAdapter47_sorted.bam PAP06_sorted.bam
  mv HTAdapter48_sorted.bam PAP07_sorted.bam
  mv HTAdapter49_sorted.bam PAP08_sorted.bam
  mv HTAdapter50_sorted.bam PAP09_sorted.bam
  mv HTAdapter51_sorted.bam PAP10_sorted.bam
  mv HTAdapter52_sorted.bam PAP11_sorted.bam
  mv HTAdapter53_sorted.bam PAP12_sorted.bam
  mv HTAdapter54_sorted.bam PAP13_sorted.bam
  mv HTAdapter55_sorted.bam PAP14_sorted.bam
  mv HTAdapter56_sorted.bam PAP15_sorted.bam
  mv HTAdapter57_sorted.bam STO01_sorted.bam
  mv HTAdapter58_sorted.bam STO02_sorted.bam
  mv HTAdapter59_sorted.bam STO03_sorted.bam
  mv HTAdapter60_sorted.bam STO04_sorted.bam
  mv HTAdapter61_sorted.bam STO05_sorted.bam
  mv HTAdapter62_sorted.bam STO06_sorted.bam
  mv HTAdapter63_sorted.bam STO07_sorted.bam
  mv HTAdapter64_sorted.bam STO08_sorted.bam
  mv HTAdapter65_sorted.bam STO09_sorted.bam
  mv HTAdapter66_sorted.bam STO10_sorted.bam
  mv HTAdapter67_sorted.bam STO11_sorted.bam
  mv HTAdapter68_sorted.bam STO12_sorted.bam
  mv HTAdapter69_sorted.bam STO13_sorted.bam
  mv HTAdapter70_sorted.bam STO14_sorted.bam
  mv HTAdapter71_sorted.bam STO15_sorted.bam
  mv HTAdapter72_sorted.bam STO17_sorted.bam
  mv HTAdapter73_sorted.bam STO18_sorted.bam
  mv HTAdapter74_sorted.bam STO19_sorted.bam
  mv HTAdapter75_sorted.bam STO20_sorted.bam
  mv HTAdapter76_sorted.bam STO21_sorted.bam
  mv HTAdapter77_sorted.bam SUM01_sorted.bam
  mv HTAdapter78_sorted.bam SUM02_sorted.bam
  mv HTAdapter79_sorted.bam SUM03_sorted.bam
  mv HTAdapter80_sorted.bam SUM04_sorted.bam
  mv HTAdapter81_sorted.bam SUM05_sorted.bam
  mv HTAdapter82_sorted.bam SUM06_sorted.bam
  mv HTAdapter83_sorted.bam SUM07_sorted.bam
  mv HTAdapter84_sorted.bam SUM08_sorted.bam
  mv HTAdapter85_sorted.bam SUM09_sorted.bam
  mv HTAdapter86_sorted.bam SUM10_sorted.bam
  mv HTAdapter87_sorted.bam SUM11_sorted.bam
  mv HTAdapter88_sorted.bam SUM12_sorted.bam
  mv HTAdapter89_sorted.bam SUM13_sorted.bam
  mv HTAdapter90_sorted.bam SUM15_sorted.bam
  mv HTAdapter91_sorted.bam SUM16_sorted.bam
  mv HTAdapter92_sorted.bam SUM17_sorted.bam
  mv HTAdapter93_sorted.bam SUM18_sorted.bam
  mv HTAdapter94_sorted.bam SUM19_sorted.bam
  mv HTAdapter95_sorted.bam SUM20_sorted.bam
  mv HTAdapter96_sorted.bam SUM21_sorted.bam
}

function readGroupsRNA {
  cd $WorkingDirectory/mappedReadsDNA
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF01_sorted.bam" O="ELF01_IDed.bam" RGPU="ELF" RGSM="ELF_54-38" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF02_sorted.bam" O="ELF02_IDed.bam" RGPU="ELF" RGSM="ELF_RA607" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF03_sorted.bam" O="ELF03_IDed.bam" RGPU="ELF" RGSM="ELF_RP54-0" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF04_sorted.bam" O="ELF04_IDed.bam" RGPU="ELF" RGSM="ELF_RP54-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF05_sorted.bam" O="ELF05_IDed.bam" RGPU="ELF" RGSM="ELF_558-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF06_sorted.bam" O="ELF06_IDed.bam" RGPU="ELF" RGSM="ELF_565-13" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF07_sorted.bam" O="ELF07_IDed.bam" RGPU="ELF" RGSM="ELF_5651-15" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF08_sorted.bam" O="ELF08_IDed.bam" RGPU="ELF" RGSM="ELF_568-3" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF09_sorted.bam" O="ELF09_IDed.bam" RGPU="ELF" RGSM="ELF_546-03" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF10_sorted.bam" O="ELF10_IDed.bam" RGPU="ELF" RGSM="ELF_547-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF11_sorted.bam" O="ELF11_IDed.bam" RGPU="ELF" RGSM="ELF_6050-3" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF12_sorted.bam" O="ELF12_IDed.bam" RGPU="ELF" RGSM="ELF_622-5" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF13_sorted.bam" O="ELF13_IDed.bam" RGPU="ELF" RGSM="ELF_636-3" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF14_sorted.bam" O="ELF14_IDed.bam" RGPU="ELF" RGSM="ELF_RA567" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF15_sorted.bam" O="ELF15_IDed.bam" RGPU="ELF" RGSM="ELF_RP54.1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="ELF16_sorted.bam" O="ELF16_IDed.bam" RGPU="ELF" RGSM="ELF_RP54" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH01_sorted.bam" O="MAH01_IDed.bam" RGPU="MAH" RGSM="MAH_35-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH02_sorted.bam" O="MAH02_IDed.bam" RGPU="MAH" RGSM="MAH_520-02" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH03_sorted.bam" O="MAH03_IDed.bam" RGPU="MAH" RGSM="MAH_539-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH04_sorted.bam" O="MAH04_IDed.bam" RGPU="MAH" RGSM="MAH_540-07" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH05_sorted.bam" O="MAH05_IDed.bam" RGPU="MAH" RGSM="MAH_541-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH06_sorted.bam" O="MAH06_IDed.bam" RGPU="MAH" RGSM="MAH_542-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH07_sorted.bam" O="MAH07_IDed.bam" RGPU="MAH" RGSM="MAH_543-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH08_sorted.bam" O="MAH08_IDed.bam" RGPU="MAH" RGSM="MAH_RA620" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH09_sorted.bam" O="MAH09_IDed.bam" RGPU="MAH" RGSM="MAH_RA641" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH10_sorted.bam" O="MAH10_IDed.bam" RGPU="MAH" RGSM="MAH_RP47.7" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH11_sorted.bam" O="MAH11_IDed.bam" RGPU="MAH" RGSM="MAH_RP47-10" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH12_sorted.bam" O="MAH12_IDed.bam" RGPU="MAH" RGSM="MAH_RP47-5" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH13_sorted.bam" O="MAH13_IDed.bam" RGPU="MAH" RGSM="MAH_RP47-8" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH14_sorted.bam" O="MAH14_IDed.bam" RGPU="MAH" RGSM="MAH_RP47-9" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAH15_sorted.bam" O="MAH15_IDed.bam" RGPU="MAH" RGSM="MAH_RP47-6" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR01_sorted.bam" O="MAR01_IDed.bam" RGPU="MAR" RGSM="MAR_RP6771" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR02_sorted.bam" O="MAR02_IDed.bam" RGPU="MAR" RGSM="MAR_RP686" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR03_sorted.bam" O="MAR03_IDed.bam" RGPU="MAR" RGSM="MAR_368- 2" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR04_sorted.bam" O="MAR04_IDed.bam" RGPU="MAR" RGSM="MAR_RA26 366-4" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR05_sorted.bam" O="MAR05_IDed.bam" RGPU="MAR" RGSM="MAR_RA379" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR06_sorted.bam" O="MAR06_IDed.bam" RGPU="MAR" RGSM="MAR_RA42 _179-4" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR07_sorted.bam" O="MAR07_IDed.bam" RGPU="MAR" RGSM="MAR_RA47_157-2T" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR08_sorted.bam" O="MAR08_IDed.bam" RGPU="MAR" RGSM="MAR_RA605" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR09_sorted.bam" O="MAR09_IDed.bam" RGPU="MAR" RGSM="MAR_RA630" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="MAR10_sorted.bam" O="MAR10_IDed.bam" RGPU="MAR" RGSM="MAR_RA88_263-3" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP01_sorted.bam" O="PAP01_IDed.bam" RGPU="PAP" RGSM="PAP_561-4" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP02_sorted.bam" O="PAP02_IDed.bam" RGPU="PAP" RGSM="PAP_562-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP03_sorted.bam" O="PAP03_IDed.bam" RGPU="PAP" RGSM="PAP_564-2" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP04_sorted.bam" O="PAP04_IDed.bam" RGPU="PAP" RGSM="PAP_569-3" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP05_sorted.bam" O="PAP05_IDed.bam" RGPU="PAP" RGSM="PAP_5010-05" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP06_sorted.bam" O="PAP06_IDed.bam" RGPU="PAP" RGSM="PAP_5021-02" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP07_sorted.bam" O="PAP07_IDed.bam" RGPU="PAP" RGSM="PAP_505-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP08_sorted.bam" O="PAP08_IDed.bam" RGPU="PAP" RGSM="PAP_516-02" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP09_sorted.bam" O="PAP09_IDed.bam" RGPU="PAP" RGSM="PAP_5241-01" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP10_sorted.bam" O="PAP10_IDed.bam" RGPU="PAP" RGSM="PAP_28-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP11_sorted.bam" O="PAP11_IDed.bam" RGPU="PAP" RGSM="PAP_RA434" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP12_sorted.bam" O="PAP12_IDed.bam" RGPU="PAP" RGSM="PAP_RA647" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP13_sorted.bam" O="PAP13_IDed.bam" RGPU="PAP" RGSM="PAP_RA648" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP14_sorted.bam" O="PAP14_IDed.bam" RGPU="PAP" RGSM="PAP_RA649" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="PAP15_sorted.bam" O="PAP15_IDed.bam" RGPU="PAP" RGSM="PAP_RP16" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO01_sorted.bam" O="STO01_IDed.bam" RGPU="STO" RGSM="STO_RP697" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO02_sorted.bam" O="STO02_IDed.bam" RGPU="STO" RGSM="STO_RP698" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO03_sorted.bam" O="STO03_IDed.bam" RGPU="STO" RGSM="STO_RA530" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO04_sorted.bam" O="STO04_IDed.bam" RGPU="STO" RGSM="STO_RA549" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO05_sorted.bam" O="STO05_IDed.bam" RGPU="STO" RGSM="STO_RP22" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO06_sorted.bam" O="STO06_IDed.bam" RGPU="STO" RGSM="STO_271" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO07_sorted.bam" O="STO07_IDed.bam" RGPU="STO" RGSM="STO_275" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO08_sorted.bam" O="STO08_IDed.bam" RGPU="STO" RGSM="STO_43-113" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO09_sorted.bam" O="STO09_IDed.bam" RGPU="STO" RGSM="STO_43-115" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO10_sorted.bam" O="STO10_IDed.bam" RGPU="STO" RGSM="STO_43-125" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO11_sorted.bam" O="STO11_IDed.bam" RGPU="STO" RGSM="STO_43-76" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO12_sorted.bam" O="STO12_IDed.bam" RGPU="STO" RGSM="STO_43-77" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO13_sorted.bam" O="STO13_IDed.bam" RGPU="STO" RGSM="STO_43-78" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO14_sorted.bam" O="STO14_IDed.bam" RGPU="STO" RGSM="STO_43-74" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO15_sorted.bam" O="STO15_IDed.bam" RGPU="STO" RGSM="STO_43-75" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO17_sorted.bam" O="STO17_IDed.bam" RGPU="STO" RGSM="STO_33-10" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO18_sorted.bam" O="STO18_IDed.bam" RGPU="STO" RGSM="STO_560-2" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO19_sorted.bam" O="STO19_IDed.bam" RGPU="STO" RGSM="STO_528-3" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO20_sorted.bam" O="STO20_IDed.bam" RGPU="STO" RGSM="STO_559-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="STO21_sorted.bam" O="STO21_IDed.bam" RGPU="STO" RGSM="STO_33-272" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM01_sorted.bam" O="SUM01_IDed.bam" RGPU="SUM" RGSM="SUM_41-103" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM02_sorted.bam" O="SUM02_IDed.bam" RGPU="SUM" RGSM="SUM_41-104" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM03_sorted.bam" O="SUM03_IDed.bam" RGPU="SUM" RGSM="SUM_RA645" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM04_sorted.bam" O="SUM04_IDed.bam" RGPU="SUM" RGSM="SUM_RP31" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM05_sorted.bam" O="SUM05_IDed.bam" RGPU="SUM" RGSM="SUM_53-255" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM06_sorted.bam" O="SUM06_IDed.bam" RGPU="SUM" RGSM="SUM_53-252" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM07_sorted.bam" O="SUM07_IDed.bam" RGPU="SUM" RGSM="SUM_53-253" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM08_sorted.bam" O="SUM08_IDed.bam" RGPU="SUM" RGSM="SUM_53-256" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM09_sorted.bam" O="SUM09_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-1" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM10_sorted.bam" O="SUM10_IDed.bam" RGPU="SUM" RGSM="SUM_644-2" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM11_sorted.bam" O="SUM11_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-13" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM12_sorted.bam" O="SUM12_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-14" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM13_sorted.bam" O="SUM13_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-5" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM15_sorted.bam" O="SUM15_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-8" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM16_sorted.bam" O="SUM16_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-9" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM17_sorted.bam" O="SUM17_IDed.bam" RGPU="SUM" RGSM="SUM_624-1" RGPL="illumina" RGLB="SeqCap2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM18_sorted.bam" O="SUM18_IDed.bam" RGPU="SUM" RGSM="SUM_625-1" RGPL="illumina" RGLB="SeqCap2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM19_sorted.bam" O="SUM19_IDed.bam" RGPU="SUM" RGSM="SUM_626-3" RGPL="illumina" RGLB="SeqCap2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM20_sorted.bam" O="SUM20_IDed.bam" RGPU="SUM" RGSM="SUM_643-2" RGPL="illumina" RGLB="SeqCap2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SUM21_sorted.bam" O="SUM21_IDed.bam" RGPU="SUM" RGSM="SUM_RP53-4" RGPL="illumina" RGLB="SeqCap2012"
}

function readGroupsRNA {
  cd $WorkingDirectory/mappedReadsRNA
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629651_sorted.bam" O="SRR629651_IDed.bam" RGPU="MAR" RGSM="MAR627-1" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629599_sorted.bam" O="SRR629599_IDed.bam" RGPU="MAH" RGSM="MAH6372" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629652_sorted.bam" O="SRR629652_IDed.bam" RGPU="MAH" RGSM="MAH6381" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629653_sorted.bam" O="SRR629653_IDed.bam" RGPU="MAR" RGSM="MAR6287" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629654_sorted.bam" O="SRR629654_IDed.bam" RGPU="MAR" RGSM="MAR276" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629655_sorted.bam" O="SRR629655_IDed.bam" RGPU="MAR" RGSM="MAR6299" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629656_sorted.bam" O="SRR629656_IDed.bam" RGPU="NAM" RGSM="NAM2064" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629657_sorted.bam" O="SRR629657_IDed.bam" RGPU="MAR" RGSM="MAR6326" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629658_sorted.bam" O="SRR629658_IDed.bam" RGPU="NAM" RGSM="NAM60603" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629659_sorted.bam" O="SRR629659_IDed.bam" RGPU="MAR" RGSM="MAR6463" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629661_sorted.bam" O="SRR629661_IDed.bam" RGPU="NAM" RGSM="NAM6193" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629660_sorted.bam" O="SRR629660_IDed.bam" RGPU="NAM" RGSM="NAM6153" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629662_sorted.bam" O="SRR629662_IDed.bam" RGPU="MAR" RGSM="MAR6099" RGPL="illumina" RGLB="RNASeq2012"
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629663_sorted.bam" O="SRR629663_IDed.bam" RGPU="NAM" RGSM="MAR6311" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629664_sorted.bam" O="SRR629664_IDed.bam" RGPU="MAR" RGSM="NAM6161" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629665_sorted.bam" O="SRR629665_IDed.bam" RGPU="MAR" RGSM="MAR6341" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629666_sorted.bam" O="SRR629666_IDed.bam" RGPU="MAH" RGSM="MAH252" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629667_sorted.bam" O="SRR629667_IDed.bam" RGPU="MAR" RGSM="MAR6503" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629668_sorted.bam" O="SRR629668_IDed.bam" RGPU="MAH" RGSM="MAH2811" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR629669_sorted.bam" O="SRR629669_IDed.bam" RGPU="MAH" RGSM="MAH6084" RGPL="illumina" RGLB="RNASeq2012" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497737_sorted.bam" O="SRR497737_IDed.bam" RGPU="ELF" RGSM="ELF523-09" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497739_sorted.bam" O="SRR497739_IDed.bam" RGPU="ELF" RGSM="ELF5171-12" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497740_sorted.bam" O="SRR497740_IDed.bam" RGPU="PAP" RGSM="PAP509-05" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497741_sorted.bam" O="SRR497741_IDed.bam" RGPU="PAP" RGSM="PAP513-03" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497742_sorted.bam" O="SRR497742_IDed.bam" RGPU="ELF" RGSM="ELF525-03" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497743_sorted.bam" O="SRR497743_IDed.bam" RGPU="ELF" RGSM="ELF544-03" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497744_sorted.bam" O="SRR497744_IDed.bam" RGPU="PAP" RGSM="PAP532-02" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497745_sorted.bam" O="SRR497745_IDed.bam" RGPU="PAP" RGSM="PAP533-07" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497746_sorted.bam" O="SRR497746_IDed.bam" RGPU="PAP" RGSM="PAP532-05" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497747_sorted.bam" O="SRR497747_IDed.bam" RGPU="ELF" RGSM="ELF5171-02" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497748_sorted.bam" O="SRR497748_IDed.bam" RGPU="ELF" RGSM="ELF524-03" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497749_sorted.bam" O="SRR497749_IDed.bam" RGPU="PAP" RGSM="PAP504-03" RGPL="illumina" RGLB="RNASeq2008" 
  java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="SRR497738_sorted.bam" O="SRR497738_IDed.bam" RGPU="ELF" RGSM="ELF523-19" RGPL="illumina" RGLB="RNASeq2008" 
}

function prepForVariantCalling {
  cd $WorkingDirectory/mappedReads"$1"
  # make .sam files list for Samtools processing
    ls | grep "_IDed.bam" |cut -d "_" -f 1 | sort | uniq  > bamList
    
    while read i;
    do
    cd $WorkingDirectory/mappedReads"$1"
    ###  remove duplicates  ###
      java -Xmx8g -jar /tools/picard-tools-2.4.1/MarkDuplicates.jar \
          INPUT=$WorkingDirectory/mappedReads"$1"/"$i"_IDed.bam \
          OUTPUT=$WorkingDirectory/GATK"$1"/"$i"_0.bam \
          METRICS_FILE=DuplicationMetrics \
          CREATE_INDEX=true \
          VALIDATION_STRINGENCY=SILENT \
          MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
          ASSUME_SORTED=TRUE \
          REMOVE_DUPLICATES=TRUE
    #index sorted bam
      samtools index $WorkingDirectory/GATK"$1"/"$i"_0.bam
    ###  Calculate Mapping Stats  ###
      samtools idxstats $WorkingDirectory/GATK"$1"/"$i"_0.bam > $WorkingDirectory/Stats"$1"/"$i"_counts.txt
      samtools flagstat $WorkingDirectory/GATK"$1"/"$i"_0.bam > $WorkingDirectory/Stats"$1"/"$i"_stats.txt
      samtools depth $WorkingDirectory/GATK"$1"/"$i"_0.bam > $WorkingDirectory/Stats"$1"/"$i"_depth.txt
  done<$WorkingDirectory/mappedReads"$1"/samList
}

function renameSortedBAMs {
# Purpose of this function: to have sampling-based naming convention be the same for both VCF and BAM files
# This aids in performing read-based masking of the fasta files

cd $WorkingDirectory/mappedReadsDNA
    mv $WorkingDirectory/mappedReadsDNA/ELF01_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRA567.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF02_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRA607.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF03_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRP540.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF04_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRP541.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF05_sorted.bam $WorkingDirectory/mappedReadsAll/ELF5581.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF06_sorted.bam $WorkingDirectory/mappedReadsAll/ELF56513.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF07_sorted.bam $WorkingDirectory/mappedReadsAll/ELF565115.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF08_sorted.bam $WorkingDirectory/mappedReadsAll/ELF5683.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF09_sorted.bam $WorkingDirectory/mappedReadsAll/ELF54603.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF10_sorted.bam $WorkingDirectory/mappedReadsAll/ELF54701.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF11_sorted.bam $WorkingDirectory/mappedReadsAll/ELF60503.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF12_sorted.bam $WorkingDirectory/mappedReadsAll/ELF6225.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF13_sorted.bam $WorkingDirectory/mappedReadsAll/ELF6363.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF14_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRA567dup.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF15_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRP54.1.bam
    mv $WorkingDirectory/mappedReadsDNA/ELF16_sorted.bam $WorkingDirectory/mappedReadsAll/ELFRP54.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH01_sorted.bam $WorkingDirectory/mappedReadsAll/MAH351.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH02_sorted.bam $WorkingDirectory/mappedReadsAll/MAH52002.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH03_sorted.bam $WorkingDirectory/mappedReadsAll/MAH53901.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH04_sorted.bam $WorkingDirectory/mappedReadsAll/MAH54007.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH05_sorted.bam $WorkingDirectory/mappedReadsAll/MAH54101.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH06_sorted.bam $WorkingDirectory/mappedReadsAll/MAH54201.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH07_sorted.bam $WorkingDirectory/mappedReadsAll/MAH54301.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH08_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRA620.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH09_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRA641.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH10_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRP47.7.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH11_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRP4710.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH12_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRP475.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH13_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRP478.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH14_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRP479.bam
    mv $WorkingDirectory/mappedReadsDNA/MAH15_sorted.bam $WorkingDirectory/mappedReadsAll/MAHRP476.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR01_sorted.bam $WorkingDirectory/mappedReadsAll/MARRP6771.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR02_sorted.bam $WorkingDirectory/mappedReadsAll/MARRP686.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR03_sorted.bam $WorkingDirectory/mappedReadsAll/MAR3682.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR04_sorted.bam $WorkingDirectory/mappedReadsAll/MAR3664.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR05_sorted.bam $WorkingDirectory/mappedReadsAll/MARRA379.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR06_sorted.bam $WorkingDirectory/mappedReadsAll/MAR1794.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR07_sorted.bam $WorkingDirectory/mappedReadsAll/MAR1572T.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR08_sorted.bam $WorkingDirectory/mappedReadsAll/MARRA605.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR09_sorted.bam $WorkingDirectory/mappedReadsAll/MARRA630.bam
    mv $WorkingDirectory/mappedReadsDNA/MAR10_sorted.bam $WorkingDirectory/mappedReadsAll/MAR2633.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP01_sorted.bam $WorkingDirectory/mappedReadsAll/PAP5614.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP02_sorted.bam $WorkingDirectory/mappedReadsAll/PAP5621.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP03_sorted.bam $WorkingDirectory/mappedReadsAll/PAP5642.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP04_sorted.bam $WorkingDirectory/mappedReadsAll/PAP5693.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP05_sorted.bam $WorkingDirectory/mappedReadsAll/PAP501005.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP06_sorted.bam $WorkingDirectory/mappedReadsAll/PAP502102.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP07_sorted.bam $WorkingDirectory/mappedReadsAll/PAP50501.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP08_sorted.bam $WorkingDirectory/mappedReadsAll/PAP51602.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP09_sorted.bam $WorkingDirectory/mappedReadsAll/PAP524101.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP10_sorted.bam $WorkingDirectory/mappedReadsAll/PAP281.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP11_sorted.bam $WorkingDirectory/mappedReadsAll/PAPRA434.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP12_sorted.bam $WorkingDirectory/mappedReadsAll/PAPRA647.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP13_sorted.bam $WorkingDirectory/mappedReadsAll/PAPRA648.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP14_sorted.bam $WorkingDirectory/mappedReadsAll/PAPRA649.bam
    mv $WorkingDirectory/mappedReadsDNA/PAP15_sorted.bam $WorkingDirectory/mappedReadsAll/PAPRP16.bam
    mv $WorkingDirectory/mappedReadsDNA/STO01_sorted.bam $WorkingDirectory/mappedReadsAll/STORP697.bam
    mv $WorkingDirectory/mappedReadsDNA/STO02_sorted.bam $WorkingDirectory/mappedReadsAll/STORP698.bam
    mv $WorkingDirectory/mappedReadsDNA/STO03_sorted.bam $WorkingDirectory/mappedReadsAll/STORA530.bam
    mv $WorkingDirectory/mappedReadsDNA/STO04_sorted.bam $WorkingDirectory/mappedReadsAll/STORA549.bam
    mv $WorkingDirectory/mappedReadsDNA/STO05_sorted.bam $WorkingDirectory/mappedReadsAll/STORP22.bam
    mv $WorkingDirectory/mappedReadsDNA/STO06_sorted.bam $WorkingDirectory/mappedReadsAll/STO271.bam
    mv $WorkingDirectory/mappedReadsDNA/STO07_sorted.bam $WorkingDirectory/mappedReadsAll/STO275.bam
    mv $WorkingDirectory/mappedReadsDNA/STO08_sorted.bam $WorkingDirectory/mappedReadsAll/STO43113.bam
    mv $WorkingDirectory/mappedReadsDNA/STO09_sorted.bam $WorkingDirectory/mappedReadsAll/STO43115.bam
    mv $WorkingDirectory/mappedReadsDNA/STO10_sorted.bam $WorkingDirectory/mappedReadsAll/STO43125.bam
    mv $WorkingDirectory/mappedReadsDNA/STO11_sorted.bam $WorkingDirectory/mappedReadsAll/STO4376.bam
    mv $WorkingDirectory/mappedReadsDNA/STO12_sorted.bam $WorkingDirectory/mappedReadsAll/STO4377.bam
    mv $WorkingDirectory/mappedReadsDNA/STO13_sorted.bam $WorkingDirectory/mappedReadsAll/STO4378.bam
    mv $WorkingDirectory/mappedReadsDNA/STO14_sorted.bam $WorkingDirectory/mappedReadsAll/STO4374.bam
    mv $WorkingDirectory/mappedReadsDNA/STO15_sorted.bam $WorkingDirectory/mappedReadsAll/STO4375.bam
    mv $WorkingDirectory/mappedReadsDNA/STO17_sorted.bam $WorkingDirectory/mappedReadsAll/STO3310.bam
    mv $WorkingDirectory/mappedReadsDNA/STO18_sorted.bam $WorkingDirectory/mappedReadsAll/STO5602.bam
    mv $WorkingDirectory/mappedReadsDNA/STO19_sorted.bam $WorkingDirectory/mappedReadsAll/STO5283.bam
    mv $WorkingDirectory/mappedReadsDNA/STO20_sorted.bam $WorkingDirectory/mappedReadsAll/STO5591.bam
    mv $WorkingDirectory/mappedReadsDNA/STO21_sorted.bam $WorkingDirectory/mappedReadsAll/STO33272.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM01_sorted.bam $WorkingDirectory/mappedReadsAll/SUM41103.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM02_sorted.bam $WorkingDirectory/mappedReadsAll/SUM41104.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM03_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRA645.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM04_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP31.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM05_sorted.bam $WorkingDirectory/mappedReadsAll/SUM53255.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM06_sorted.bam $WorkingDirectory/mappedReadsAll/SUM53252.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM07_sorted.bam $WorkingDirectory/mappedReadsAll/SUM53253.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM08_sorted.bam $WorkingDirectory/mappedReadsAll/SUM53256.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM09_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP531.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM10_sorted.bam $WorkingDirectory/mappedReadsAll/SUM6442.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM11_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP5313.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM12_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP5314.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM13_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP535.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM15_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP538.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM16_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP539.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM17_sorted.bam $WorkingDirectory/mappedReadsAll/SUM6241.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM18_sorted.bam $WorkingDirectory/mappedReadsAll/SUM6251.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM19_sorted.bam $WorkingDirectory/mappedReadsAll/SUM6263.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM20_sorted.bam $WorkingDirectory/mappedReadsAll/SUM6432.bam
    mv $WorkingDirectory/mappedReadsDNA/SUM21_sorted.bam $WorkingDirectory/mappedReadsAll/SUMRP534.bam

cd $WorkingDirectory/mappedReadsRNA
    mv $WorkingDirectory/mappedReadsRNA/SRR497737_sorted.bam $WorkingDirectory/mappedReadsAll/ELF52309.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497738_sorted.bam $WorkingDirectory/mappedReadsAll/ELF52319.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497739_sorted.bam $WorkingDirectory/mappedReadsAll/ELF517112.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497740_sorted.bam $WorkingDirectory/mappedReadsAll/PAP50905.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497741_sorted.bam $WorkingDirectory/mappedReadsAll/PAP51303.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497742_sorted.bam $WorkingDirectory/mappedReadsAll/ELF52503.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497743_sorted.bam $WorkingDirectory/mappedReadsAll/ELF54403.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497744_sorted.bam $WorkingDirectory/mappedReadsAll/PAP53202.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497745_sorted.bam $WorkingDirectory/mappedReadsAll/PAP53307.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497746_sorted.bam $WorkingDirectory/mappedReadsAll/PAP53205.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497747_sorted.bam $WorkingDirectory/mappedReadsAll/ELF517102.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497748_sorted.bam $WorkingDirectory/mappedReadsAll/ELF52403.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR497749_sorted.bam $WorkingDirectory/mappedReadsAll/PAP50403.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629599_sorted.bam $WorkingDirectory/mappedReadsAll/MAH6372.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629651_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6271.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629652_sorted.bam $WorkingDirectory/mappedReadsAll/MAH6381.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629653_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6287.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629654_sorted.bam $WorkingDirectory/mappedReadsAll/MAR276.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629655_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6299.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629656_sorted.bam $WorkingDirectory/mappedReadsAll/NAM2064.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629657_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6326.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629658_sorted.bam $WorkingDirectory/mappedReadsAll/NAM60603.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629659_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6463.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629660_sorted.bam $WorkingDirectory/mappedReadsAll/NAM6153.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629661_sorted.bam $WorkingDirectory/mappedReadsAll/NAM6193.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629662_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6099.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629663_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6311.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629664_sorted.bam $WorkingDirectory/mappedReadsAll/NAM6161.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629665_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6341.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629666_sorted.bam $WorkingDirectory/mappedReadsAll/MAH252.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629667_sorted.bam $WorkingDirectory/mappedReadsAll/MAR6503.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629668_sorted.bam $WorkingDirectory/mappedReadsAll/MAH2811.bam
    mv $WorkingDirectory/mappedReadsRNA/SRR629669_sorted.bam $WorkingDirectory/mappedReadsAll/MAH6084.bam

}

# -- PREPARE FOR SNP CALLING -- ##
# •••••••••••• Functions for obtaining SNP convergence in bqsr ••••••••••••• #
# •••••••••• The input parameter is the current replicate number ••••••••••• #
# ••••••• This approach is called "bootstrapping" by GATK developers ••••••• #
function indexReference {
  # INDEX REF TO USE FOR SNP CALLING
    samtools faidx $WorkingDirectory/References/TelagGenome.fasta
    java -Xmx8g -jar /tools/picard-tools-2.4.1/CreateSequenceDictionary.jar \
          R= $WorkingDirectory/References/TelagGenome.fasta \
          O= $WorkingDirectory/References/TelagGenome.dict 
}

function use-HaplotypeCaller {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" HaplotypeCaller \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -I "$i"_"$1".bam \
    -O "$i"_"$1".g.vcf \
    -ERC GVCF
  echo "$i"$'\t'"$i""_""$1"".g.vcf" >> cohort.sample_map_"$1"
done<$WorkingDirectory/mappedReads"$2"/samList2
}

function get-just-SNPs {
    # Combine GVCFs
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" GenomicsDBImport \
    	--sample-name-map cohort.sample_map_"$1" \
    	--genomicsdb-workspace-path SNP_database_"$1" \
    	--reader-threads 5 \
    	--intervals genome.intervals
    # Joint genotyping
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" GenotypeGVCFs \
    	-R $WorkingDirectory/References/TelagGenome.fasta \
    	-V gendb://SNP_database_"$1" \
    	-O genotyped_"$1".vcf 
    # Get SNPs
    /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" SelectVariants \
	-R $WorkingDirectory/References/TelagGenome.fasta \
	-V genotyped_"$1".vcf \
	--select-type-to-include SNP \
	-O JustSNPs_"$1".vcf
}

function use-BaseRecalibrator {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" BaseRecalibrator \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -I "$i"_"$1".bam \
    --known-sites filtered_0.vcf \
    -O "$i"_recalibration_"$1".table
done<$WorkingDirectory/mappedReads"$2"/samList
}

function use-AnalyzeCovariates {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" AnalyzeCovariates \
    -before "$i"_recalibration_"$1".table \
    -after "$i"_recalibration_"$2".table \
    -plots "$i"_recalComparison_"$1".pdf \
    -csv "$i"_recalComparison_"$1".csv
done<$WorkingDirectory/mappedReads"$3"/samList
}

function use-BQSR {
while read i
do
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" ApplyBQSR \
    -R $WorkingDirectory/References/TelagGenome.fasta \
    -I "$i"_"$1".bam \
    --bqsr-recal-file "$i"_recalibration_"$1".table \
    -O "$i"_"$2".bam
done<$WorkingDirectory/mappedReads"$3"/samList
}

function combine-VCF {
    cp $WorkingDirectory/GATKDNA/JustSNPs_1.vcf $WorkingDirectory/variantFiltration/JustSNPs_DNA.vcf
    cp $WorkingDirectory/GATKRNA/JustSNPs_1.vcf $WorkingDirectory/variantFiltration/JustSNPs_RNA.vcf
    cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/variantFiltration
    # Change names:
    # DNA
    sed -i.bak "s/SRR497737/ELF52309/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497738/ELF52319/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497739/ELF517112/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497740/PAP50905/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497741/PAP51303/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497742/ELF52503/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497743/ELF54403/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497744/PAP53202/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497745/PAP53307/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497746/PAP53205/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497747/ELF517102/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497748/ELF52403/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR497749/PAP50403/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629599/MAH6372/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629651/MAR6271/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629652/MAR6287/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629653/MAR276/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629654/MAR6299/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629655/NAM2064/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629656/MAR6326/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629657/MAR6326dup/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629658/NAM60603/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629659/MAR6463/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629660/NAM6193/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629661/NAM6153/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629662/MAR6099/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629663/NAM6161/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629664/MAR6311/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629665/MAR6341/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629666/MAH252/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629667/MAR6503/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629668/MAH2811/" JustSNPs_RNA.vcf
    sed -i.bak "s/SRR629669/MAH6084/" JustSNPs_RNA.vcf

    # RNA
    sed -i.bak "s/ELF01/ELFRA567/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF02/ELFRA607/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF03/ELFRP540/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF04/ELFRP541/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF05/ELF5581/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF06/ELF56513/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF07/ELF565115/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF08/ELF5683/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF09/ELF54603/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF10/ELF54701/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF11/ELF60503/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF12/ELF6225/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF13/ELF6363/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF14/ELFRA567dup/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF15/ELFRP54.1/" JustSNPs_DNA.vcf
    sed -i.bak "s/ELF16/ELFRP54/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH01/MAH351/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH02/MAH52002/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH03/MAH53901/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH04/MAH54007/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH05/MAH54101/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH06/MAH54201/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH07/MAH54301/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH08/MAHRA620/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH09/MAHRA641/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH10/MAHRP47.7/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH11/MAHRP4710/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH12/MAHRP475/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH13/MAHRP478/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH14/MAHRP479/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAH15/MAHRP476/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR01/MARRP6771/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR02/MARRP686/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR03/MAR3682/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR04/MAR3664/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR05/MARRA379/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR06/MAR1794/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR07/MAR1572T/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR08/MARRA605/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR09/MARRA630/" JustSNPs_DNA.vcf
    sed -i.bak "s/MAR10/MAR2633/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP01/PAP5614/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP02/PAP5621/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP03/PAP5642/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP04/PAP5693/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP05/PAP501005/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP06/PAP502102/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP07/PAP50501/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP08/PAP51602/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP09/PAP524101/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP10/PAP281/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP11/PAPRA434/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP12/PAPRA647/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP13/PAPRA648/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP14/PAPRA649/" JustSNPs_DNA.vcf
    sed -i.bak "s/PAP15/PAPRP16/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO01/STORP697/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO02/STORP698/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO03/STORA530/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO04/STORA549/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO05/STORP22/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO06/STO271/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO07/STO275/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO08/STO43113/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO09/STO43115/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO10/STO43125/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO11/STO4376/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO12/STO4377/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO13/STO4378/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO14/STO4374/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO15/STO4375/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO17/STO3310/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO18/STO5602/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO19/STO5283/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO20/STO5591/" JustSNPs_DNA.vcf
    sed -i.bak "s/STO21/STO33272/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM01/SUM41103/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM02/SUM41104/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM03/SUMRA645/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM04/SUMRP31/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM05/SUM53255/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM06/SUM53252/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM07/SUM53253/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM08/SUM53256/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM09/SUMRP531/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM10/SUM6442/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM11/SUMRP5313/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM12/SUMRP5314/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM13/SUMRP535/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM15/SUMRP538/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM16/SUMRP539/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM17/SUM6241/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM18/SUM6251/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM19/SUM6263/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM20/SUM6432/" JustSNPs_DNA.vcf
    sed -i.bak "s/SUM21/SUMRP534/" JustSNPs_DNA.vcf
    bgzip JustSNPs_DNA.vcf
    bgzip JustSNPs_RNA.vcf
    bcftools index JustSNPs_DNA.vcf.gz
    bcftools index JustSNPs_RNA.vcf.gz
    bcftools merge JustSNPs_RNA.vcf.gz JustSNPs_DNA.vcf.gz -O v -o Merged.vcf
}

function removeRNAedits {
  vcftools --vcf "$1".vcf --non-ref-af 1 \
  --recode --recode-INFO-all --out Substitutions.vcf
  mv Substitutions.vcf.recode.vcf Substitutions.vcf
  grep -v "#" Substitutions.vcf | cut -f 1,2 > Substitutions.txt
  vcftools --vcf "$1".vcf --exclude-positions Substitutions.txt \
  --recode --recode-INFO-all --out removedRNAedits
  mv removedRNAedits.recode.vcf removedRNAedits.vcf
  echo "RNA edits removed variants: $(grep -v "^#" removedRNAedits.vcf | wc -l)" >> Log.txt
}

function getNetworkFasta {
 # Extract genes from exons used for probe design
 cd $WorkingDirectory/References
 python $pythonScripts/filterExons.py Exons_2021.fa "$1".txt "$1"TargetGenes.fa Extract"$1"GenesFromExons_log.txt
}


function probes2gff {
 ## ANNOTATE FILTERED VARIANT FILES FOR SEQCAP DATA
 cd $WorkingDirectory/References
 #+ COMPLETED # Make blast database from T. elegans genome
 #+ COMPLETED makeblastdb -in TelagGenome.fasta -parse_seqids -dbtype nucl -out Genome.db
 #+ COMPLETED # Use Blast with Exons_2021.fa (the exons used for probe design) to filter the genome
 #+ COMPLETED blastn -db Genome.db -query "$1" -outfmt "7 qseqid sseqid evalue qstart qend sstart send" -out BlastResults_"$2".txt
 #+ COMPLETED # Delete "^#" lines from blast output
 #+ COMPLETED cp BlastResults_"$2".txt BlastResults_"$2"_original.txt
 #+ COMPLETED sed -i.bak '/^#/d' BlastResults_"$2".txt
 #+ COMPLETED sed -i.bak "s/ref|//" BlastResults_"$2".txt
 #+ COMPLETED sed -i.bak "s/|//" BlastResults_"$2".txt
 # Use filtered genome results (blast output) to pull out targeted genes and create filtered gff
 python $pythonScripts/shrinkGFF.py BlastResults_"$2".txt TelagGenome.gff "$2"_CapturedGenes.gff "$2"_CapturedExons.gff "$2"_CapturedCDS.gff Pull"$2"CapturedGenes_log.txt
 # Use bedops to convert gff to bed
 gff2bed < "$2"_CapturedGenes.gff > "$2"_CapturedGenes.bed
 gff2bed < "$2"_CapturedExons.gff > "$2"_CapturedExons.bed
 gff2bed < "$2"_CapturedCDS.gff > "$2"_CapturedCDS.bed
 # bgzip bed file
 bgzip -f "$2"_CapturedGenes.bed "$2"_CapturedGenes.bed.gz
 bgzip -f "$2"_CapturedExons.bed "$2"_CapturedExons.bed.gz
 bgzip -f "$2"_CapturedCDS.bed "$2"_CapturedCDS.bed.gz
 # tabix index .bed.gz file
 tabix -f -p bed "$2"_CapturedGenes.bed.gz
 tabix -f -p bed "$2"_CapturedExons.bed.gz
 tabix -f -p bed "$2"_CapturedCDS.bed.gz
}

function annotateVariants {
  # Annotate SNP file make sure "$1".vcf is in this directory
  cd $WorkingDirectory/variantFiltration
  bcftools annotate \
  	-a $WorkingDirectory/References/"$2"_CapturedGenes.bed.gz \
  	-c CHROM,FROM,TO,GENE \
        -o "$2"_Annotated_Init.vcf \
  	-O v \
  	-h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  	"$1".vcf
  awk '/^#|GENE=/' "$2"_Annotated_Init.vcf > "$2"_Annotated.vcf
  echo "Annotated "$2" variants: $(grep -v "^#" "$2"_Annotated.vcf | wc -l)" >> AnnotationLog.txt
}  

function getSpecificVariants {
  # Annotate SNP file make sure "$1".vcf is in this directory
  cd $WorkingDirectory/variantFiltration
  bcftools annotate -x INFO/GENE "$1"_HardFilterStep3.vcf > "$1"_Deannotated.vcf
  bcftools annotate \
  	-a $WorkingDirectory/References/"$1"_Captured"$2".bed.gz \
  	-c CHROM,FROM,TO,GENE \
        -o "$1"_"$2"_Init.vcf \
  	-O v \
  	-h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
  	"$1"_Deannotated.vcf
  awk '/^#|GENE=/' "$1"_"$2"_Init.vcf > "$1"_"$2".vcf
  echo "CDS "$1" variants: $(grep -v "^#" "$1"_"$2".vcf | wc -l)" >> Log.txt
}

function filterByPopulation {
  cd $WorkingDirectory/variantFiltration
  #Create list for each population
  echo "ELF" >> Pops; echo "MAH" >> Pops; echo "MAR" >> Pops; echo "NAM" >> Pops; echo "PAP" >> Pops; echo "STO" >> Pops; echo "SUM" >> Pops
  while read i
  do
    #Create VCF for each population  
    bcftools view --samples-file $WorkingDirectory/References/"$i".txt "$1" > "$2"_"$i".vcf
    #Get number of individuals for each vcf
    n="$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' "$2"_"$i".vcf)"
    min=$(expr $n - 1)
    #Filter missing data- 
    vcftools --max-missing-count $min --vcf "$2"_"$i".vcf  --recode --recode-INFO-all --out "$2"_"$i"_popFiltered.vcf
    mv "$2"_"$i"_popFiltered.vcf.recode.vcf "$2"_"$i"_popFiltered.vcf
    bgzip "$2"_"$i"_popFiltered.vcf
    bcftools index -f "$2"_"$i"_popFiltered.vcf.gz
  done<Pops
  bcftools merge "$2"*popFiltered.vcf.gz -O v -o "$2"_popFiltered.vcf
  echo "Filtered by Population "$2" variants: $(grep -v "^#" "$2"_popFiltered.vcf | wc -l)" >> Log.txt
}

function plotVariants {
source /home/rlk0015/miniconda3/etc/profile.d/conda.sh
conda activate vcfEnv
python $pythonScripts/vcf2table.py "$1" "$1"_table
~/SeqCap/RScripts/plotvcftable.R -I "$1"_table -O "$1"_plot.pdf
}

function initial-VariantFiltration {
# Filter variants
	# Provide expression for filtration:
	  # FS = FisherStrand *NOT COMMON IN GATK4.1.7.0
	    # Phred-scaled probability strand bias exists at site 
	  # QD = quality of depth *NOT COMMON IN GATK4.1.7.0
	    # Variant confidence divided by raw depth 
	  # DP = depth of coverage *NOT COMMON IN GATK4.1.7.0
	    # Raw depth. I don't use here, because QD accounts for this 
	  # MQ = RMS mapping quality *NOT COMMON IN GATK4.1.7.0
	    # Root mean square mq over all the reads at the site 
	  # MQRankSum = Mapping quality rank sum test *NOT COMMON IN GATK4.1.7.0
	    # Compares mapping qualities of reads supporting the ref allele to those supporting the alt allele. 
	  # ReadPosRankSum = Rank sum test used to assess site position *NOT COMMON IN GATK4.1.7.0
	    # It compares whether the positions of the reference and alternate alleles are different within the reads 
	  # AB = Allele Balance
	    # Depth of covereage for each allele per sample generalized over all samples
	  # MQ0 = Map Quality 0
	    # Number of reads with map quality 0
	  # SOR = StrandOddsRatio
	    # Estimate strand bias without penalizing reads that occur at the end of exons *FS is biased to penalize
  /tools/gatk-4.1.7.0/gatk --java-options "-Xmx16g" VariantFiltration \
	-R $WorkingDirectory/References/TelagGenome.fasta \
	-V "$1" \
        -O filtered_"$2"Init.vcf \
	--filter-name "SOR" \
	--filter-expression "SOR > 3.0" \
	--filter-name "QD" \
	--filter-expression "QD < 2.0" \
	--filter-name "MQ" \
	--filter-expression "MQ < 40.0" \
	--filter-name "MQRankSum" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "FS" \
	--filter-expression "FS > 60.0" \
	--filter-name "ReadPosRankSum" \
	--filter-expression "ReadPosRankSum < -5.0"
  # don't do the following for filtered_0.vcf, just get rid of the Init suffix
  awk '/^#/||$7=="PASS"' filtered_"$2"Init.vcf > filtered_"$2".vcf
  echo "Initial filtration "$2" variants: $(grep -v "^#" filtered_"$2".vcf | wc -l)" >> Log.txt
} 

function hard-VariantFiltration {
  # Step 0: Get rid of unwanted individuals (T. sirtalis, duplicates, and siblings)
    vcftools --remove Full_IndividualsToRemove --vcf filtered_"$1".vcf --recode --out "$2"_HardFilterStep0.vcf
    mv "$2"_HardFilterStep0.vcf.recode.vcf "$2"_HardFilterStep0.vcf
  # Step 1: Get rid of low-quality (mean) genotyping:
    bcftools view  -i  'MIN(FMT/GQ>20)' "$2"_HardFilterStep0.vcf > "$2"_HardFilterStep1.vcf
    echo "Genotype Quality filtration $2 variants: $(grep -v "^#" "$2"_HardFilterStep1.vcf | wc -l)" >> Log.txt
  # Step 2: Get rid of multiallelic SNPs (more than 2 alleles):
    bcftools view -m2 -M2 -v snps "$2"_HardFilterStep1.vcf > "$2"_HardFilterStep2.vcf
    echo "Multiallelic filtration $2 variants: $(grep -v "^#" "$2"_HardFilterStep2.vcf | wc -l)" >> Log.txt
  # Step 3: Filter individuals by quality
    vcftools --minGQ 20 --vcf "$2"_HardFilterStep2.vcf --recode --recode-INFO-all --out "$2"_HardFilterStep3.vcf
    mv "$2"_HardFilterStep3.vcf.recode.vcf "$2"_HardFilterStep3.vcf
    echo "Indiv GQ filtration $2 variants:  $(grep -v "^#" "$2"_HardFilterStep3.vcf | wc -l)" >> Log.txt
  # Step 4: Get rid of low-depth individuals per site
    bcftools view  -i  'MIN(FMT/DP>9)' "$2"_HardFilterStep3.vcf > "$2"_HardFilterStep4.vcf 
    echo "Genotype Quality filtration $2 variants: $(grep -v "^#" "$2"_HardFilterStep4.vcf | wc -l)" >> Log.txt
  # Step 5: Get rid of low-frequency alleles- here just singletons:
    vcftools --mac 2 --vcf "$2"_HardFilterStep4.vcf --recode --recode-INFO-all --out "$2"_HardFilterStep5.vcf
    mv "$2"_HardFilterStep5.vcf.recode.vcf "$2"_HardFilterStep5.vcf
    echo "Singleton filtration $2 variants: $(grep -v "^#" "$2"_HardFilterStep5.vcf | wc -l)" >> Log.txt

}

### -- If needed, here is code for removing individuals that would need to be adapted
# To remove bad individuals:
## vcftools --vcf yourvcfhere.vcf --depth --out Prefix
## vcftools --vcf yourvcf.vcf --remove removetheseind.txt --recode --recode-INFO-all --out gooddepth
## -- Examine Unfiltered Variants

