#/bin/sh

#################
# Purpose: Process 2008 Garter Snake Heat Stress RNAseq data 
  # Downoload from NCBI SRA
  # Determine Quality with Fastqc
  # Clean with Trimmomatic
  # Map to reference with HiSat2
  # count transcripts with Stringtie
# Author: Tonia Schwartz
#################

#-- Auburn University High Performance and Parallel Computing
#-- Hopper Cluster Sample Job Submission Script

#-- This script provides the basic scheduler directives you
#-- can use to submit a job to the Hopper scheduler.
#-- Other than the last two lines it can be used as-is to
#-- send a single node job to the cluster. Normally, you
#-- will want to modify the #PBS directives below to reflect
#-- your workflow...

####-- For convenience, give your job a name

#PBS -N RNAseq2008_ToGeneome

#-- Provide an estimated wall time in which to run your job
#-- The format is DD:HH:MM:SS.  


#PBS -l walltime=01:00:00:00 

#-- Indicate if\when you want to receive email about your job
#-- The directive below sends email if the job is (a) aborted, 
#-- when it (b) begins, and when it (e) ends

#PBS -m abe tss0019@auburn.edu

#-- Inidicate the working directory path to be used for the job.
#-- If the -d option is not specified, the default working directory 
#-- is the home directory. Here, we set the working directory
#-- current directory

#PBS -d /home/tss0019/GarterSnakes

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional

#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.

#PBS -l nodes=1:ppn=20

#-- Now issue the commands that you want to run on the compute nodes.

#-- With the -V option, you can load any software modules
#-- either before submitting, or in the job submission script.

#-- You should modify the lines below to reflect your own
#-- workflow...

#module load <myprogram_modulefile>

module load gnu_parallel/201612222
module load fastqc/11.5
module load trimmomatic/0.37

module load hisat/2.1.0
#module load stringtie/1.3.2d
module load stringtie/1.3.3b

#module load python/2.7.11
module load python/2.7.14

module load gcc/5.1.0
module load samtools/1.6
module load bcftools/1.3.1
module load gffread/2
module load gffcompare/1
#./myprogram <parameters>

#--  After saving this script, you can submit your job to the queue
#--  with...

#--  qsub sample_job.sh
##########################################


#PBS -j oe
#PBS -q debug


##  Set the stack size to unlimited
ulimit -s unlimited

# #Turn echo on so all commands are echoed in the output log
set -x

## Make directories
mkdir /scratch/tss0019/Telegans/Mapping/2008Data
mkdir /scratch/schwartz_T_lab/GarterSnake/RNAseq_2008Data/RawData
mkdir /scratch/schwartz_T_lab/GarterSnake/RNAseq_2008Data/CleanedData
mkdir /scratch/schwartz_T_lab/GarterSnake/RNAseq_2008Data/FastQC_Quality
mkdir /scratch/tss0019/Telegans/Mapping/2008Data/Results


## Define Directory Variables
DATADIR=/scratch/schwartz_T_lab/GarterSnake/RNAseq_2008Data/RawData
DATADIR_CLEAN=/scratch/schwartz_T_lab/GarterSnake/RNAseq_2008Data/CleanedData
DATADIR_QUALITY=/scratch/schwartz_T_lab/GarterSnake/RNAseq_2008Data/FastQC_Quality
OUTDIR=/scratch/tss0019/Telegans/Mapping/2008Data/Results

## This is the location of the Reference Genome and annotation
REFDIR=/scratch/tss0019/Telegans/Mapping/Genome2019_rThaEle1
REF=TelagGenome

###########  Download the RNAseq data

cd $DATADIR
fastq-dump -F --split-files SRR497737
fastq-dump -F --split-files SRR497738
fastq-dump -F --split-files SRR497739
fastq-dump -F --split-files SRR497740
fastq-dump -F --split-files SRR497741
fastq-dump -F --split-files SRR497742
fastq-dump -F --split-files SRR497743
fastq-dump -F --split-files SRR497744
fastq-dump -F --split-files SRR497745
fastq-dump -F --split-files SRR497746
fastq-dump -F --split-files SRR497747
fastq-dump -F --split-files SRR497748
fastq-dump -F --split-files SRR497749

######### Assess data quality

cd $DATADIR

##  Run fastqc on the All files
fastqc *.fastq
mv *html  $DATADIR_QUALITY
mv *zip   $DATADIR_QUALITY


############  Clean the reads
############ Trimmomatic #############
############  Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
#MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
#SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across
#requiredQuality: specifies the average quality required.
## -threads  is the option to define the number of threads (cores) to use. 
        #For this to be effective you need to request those cores at submission

### Make a while loop to work through a list of the files
  # First, create list of fastq files
  # Example file name: SRR497737_1.fastq
  # grab all fasta files, cut on the underscore, use only the first of the cuts, sort, use unique
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq  > list   #should list SRR497737 (for example)

while read i;
do
  ## count number of lines before cleaning
wc -l "$i"_1.fastq  >>  LineCount_Raw.txt
  ## trim 
java -jar /tools/trimmomatic-0.37/bin/trimmomatic.jar SE -threads 6 -phred33  "$i"_1.fastq  $DATADIR_CLEAN/"$i"_cleaned.fastq  LEADING:20 TRAILING:20 SLIDINGWINDOW:6:20 MINLEN:36

done<list


###############  Prepare the reference Index for HiSat2
	## Move to $REFDIR
cd $REFDIR
  ## convert the annotation file from .gff to .gtf
gffread "$REFDIR"/"$REF".gff -T -o "$REFDIR"/"$REF".gtf

extract_splice_sites.py "$REFDIR"/"$REF".gtf > "$REFDIR"/"$REF".ss
extract_exons.py "$REFDIR"/"$REF".gtf > "$REFDIR"/"$REF".exon

	###Create a HISAT2 index
hisat2-build --ss "$REFDIR"/"$REF".ss --exon "$REFDIR"/"$REF".exon "$REFDIR"/"$REF".fasta "$REF"_index


##############  Map and Count the 2008 SE data
# Move to the 2008 SE data directory
cd $DATADIR_CLEAN
wc -l *.fastq  >>  LineCount_Cleaned.txt

  ## Create list of fastq files to map
  # Example: SRR497747_cleaned.fastq
  # grab all cleaned fasta files, cut on the underscore, use only the first of the cuts, sort, use unique
ls | grep ".fastq" |cut -d "_" -f 1| sort | uniq > list   

## Move to the Output directory 
cd $OUTDIR

  # copy the list of unique ids from the cleaned files to map
  # This is from the 2008 SE directory
cp $DATADIR_CLEAN/list . 

mkdir ballgown

while read i;
do
     ## Map, HiSat  -p indicates 20 processors, --dta reports alignments for StringTieq
hisat2 -p 20 --dta -x $REFDIR/"$REF"_index -U $DATADIR_CLEAN/"$i"_cleaned.fastq -S "$i".sam

    ### view: convert the SAM file into a BAM file  -bS: BAM is the binary format corresponding to the SAM text format.
    ### sort: convert the BAM file to a sorted BAM file.
    ### Example Input: HS06_GATCAG_All.sam; Output: HS06_GATCAG_sorted.bam
samtools view -@ 20 -bS "$i".sam  | samtools sort -@ 20 -o  "$i"_sorted.bam   
    ### make directory for sample to hold the transcript counts
mkdir ballgown/"$i"
  # Original: This will make transcripts using the reference geneome as a guide for each sorted.bam
  # eBG: This will run stringtie once and ONLY use the Ref annotation (not make new transcripts)
stringtie -p 20 -e -B -G $REFDIR/"$REF".gtf -o ballgown/"$i"/"$i".gtf  -l "$i"  $OUTDIR/"$i"_sorted.bam

done<list

###################  Prepare the Data for Differential Gene Expression Analysis with EdgeR/Limma
### The PrepDE.py is a python script that converts the files in your ballgown folder to a count matrix in .csv format
# move out a directory so this script can work on the output directory that contains all the Ballgown subdirectories
cd ..
python /home/scripts/PrepDE.py $OUTDIR
    ##  The final count matrix is a single .csv file containing counts for all individuals and all annotated genes.
    ##  Use this file for DGE