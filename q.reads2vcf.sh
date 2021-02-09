#!/bin/sh

#Give job a name
#PBS -N PE_Genome

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=2:ppn=10,walltime=40:00:00 
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

module load fastqc/11.5
module load gnu-parallel/20160322 
module load trimmomatic/0.36
module load bwa/0.7.15
module load samtools/1.3.1
module load picard/2.4.1
module load python/3.5.2

# #create working directory
# mkdir -p /scratch/rlk0015/Telag/Dec2017Results/PE/Genome_WorkingDirectory
# echo "RLK_report: directory created: /scratch/rlk0015/Telag/Dec2017Results/PE/Genome_WorkingDirectory"
# 
# # copy over clean PE fastq files
# cd /scratch/rlk0015/Telag/Dec2017Results/PE/PE_CleanReads
# ls *_paired.fastq | parallel -j+0  --eta 'cp {} /scratch/rlk0015/Telag/Dec2017Results/PE/Genome_WorkingDirectory'
# echo "RLK_report: clean PE fastq files copied to /scratch/rlk0015/Telag/Dec2017Results/PE/Genome_WorkingDirectory"

# move to working directory
cd /scratch/rlk0015/Telag/Dec2017Results/PE/Genome_WorkingDirectory

# #copy the reference to your current directory
# cp /home/rlk0015/SeqCap/code/References/Genome.fa .
echo "RLK_report: reference copy and unzip complete"
# index reference
bwa index -p genome -a bwtsw Genome.fa
echo "RLK_report: reference index complete"

# create list with each paired individual
ls | grep "_paired.fastq" | cut -d "_" -f 1 | sort | uniq > pairedMapList

# while loop through the names in pairedMapList
while read i
do
# map to ref
bwa mem -t 4 -M genome "$i"*_R1_paired.fastq "$i"*_R2_paired.fastq > "$i"_doubles.sam
done<pairedMapList
echo "RLK_report: map complete"

# map: Tonia's method
mv HTAdapter1_doubles.sam ELF01_doubles.sam
mv HTAdapter2_doubles.sam ELF02_doubles.sam
mv HTAdapter3_doubles.sam ELF03_doubles.sam
mv HTAdapter4_doubles.sam ELF04_doubles.sam
mv HTAdapter5_doubles.sam ELF05_doubles.sam
mv HTAdapter6_doubles.sam ELF06_doubles.sam
mv HTAdapter7_doubles.sam ELF07_doubles.sam
mv HTAdapter8_doubles.sam ELF08_doubles.sam
mv HTAdapter9_doubles.sam ELF09_doubles.sam
mv HTAdapter10_doubles.sam ELF10_doubles.sam
mv HTAdapter11_doubles.sam ELF11_doubles.sam
mv HTAdapter12_doubles.sam ELF12_doubles.sam
mv HTAdapter13_doubles.sam ELF13_doubles.sam
mv HTAdapter14_doubles.sam ELF14_doubles.sam
mv HTAdapter15_doubles.sam ELF15_doubles.sam
mv HTAdapter16_doubles.sam ELF16_doubles.sam
mv HTAdapter17_doubles.sam MAH01_doubles.sam
mv HTAdapter18_doubles.sam MAH02_doubles.sam
mv HTAdapter19_doubles.sam MAH03_doubles.sam
mv HTAdapter20_doubles.sam MAH04_doubles.sam
mv HTAdapter21_doubles.sam MAH05_doubles.sam
mv HTAdapter22_doubles.sam MAH06_doubles.sam
mv HTAdapter23_doubles.sam MAH07_doubles.sam
mv HTAdapter24_doubles.sam MAH08_doubles.sam
mv HTAdapter25_doubles.sam MAH09_doubles.sam
mv HTAdapter26_doubles.sam MAH10_doubles.sam
mv HTAdapter27_doubles.sam MAH11_doubles.sam
mv HTAdapter28_doubles.sam MAH12_doubles.sam
mv HTAdapter29_doubles.sam MAH13_doubles.sam
mv HTAdapter30_doubles.sam MAH14_doubles.sam
mv HTAdapter31_doubles.sam MAH15_doubles.sam
mv HTAdapter32_doubles.sam MAR01_doubles.sam
mv HTAdapter33_doubles.sam MAR02_doubles.sam
mv HTAdapter34_doubles.sam MAR03_doubles.sam
mv HTAdapter35_doubles.sam MAR04_doubles.sam
mv HTAdapter36_doubles.sam MAR05_doubles.sam
mv HTAdapter37_doubles.sam MAR06_doubles.sam
mv HTAdapter38_doubles.sam MAR07_doubles.sam
mv HTAdapter39_doubles.sam MAR08_doubles.sam
mv HTAdapter40_doubles.sam MAR09_doubles.sam
mv HTAdapter41_doubles.sam MAR10_doubles.sam
mv HTAdapter42_doubles.sam PAP01_doubles.sam
mv HTAdapter43_doubles.sam PAP02_doubles.sam
mv HTAdapter44_doubles.sam PAP03_doubles.sam
mv HTAdapter45_doubles.sam PAP04_doubles.sam
mv HTAdapter46_doubles.sam PAP05_doubles.sam
mv HTAdapter47_doubles.sam PAP06_doubles.sam
mv HTAdapter48_doubles.sam PAP07_doubles.sam
mv HTAdapter49_doubles.sam PAP08_doubles.sam
mv HTAdapter50_doubles.sam PAP09_doubles.sam
mv HTAdapter51_doubles.sam PAP10_doubles.sam
mv HTAdapter52_doubles.sam PAP11_doubles.sam
mv HTAdapter53_doubles.sam PAP12_doubles.sam
mv HTAdapter54_doubles.sam PAP13_doubles.sam
mv HTAdapter55_doubles.sam PAP14_doubles.sam
mv HTAdapter56_doubles.sam PAP15_doubles.sam
mv HTAdapter57_doubles.sam STO01_doubles.sam
mv HTAdapter58_doubles.sam STO02_doubles.sam
mv HTAdapter59_doubles.sam STO03_doubles.sam
mv HTAdapter60_doubles.sam STO04_doubles.sam
mv HTAdapter61_doubles.sam STO05_doubles.sam
mv HTAdapter62_doubles.sam STO06_doubles.sam
mv HTAdapter63_doubles.sam STO07_doubles.sam
mv HTAdapter64_doubles.sam STO08_doubles.sam
mv HTAdapter65_doubles.sam STO09_doubles.sam
mv HTAdapter66_doubles.sam STO10_doubles.sam
mv HTAdapter67_doubles.sam STO11_doubles.sam
mv HTAdapter68_doubles.sam STO12_doubles.sam
mv HTAdapter69_doubles.sam STO13_doubles.sam
mv HTAdapter70_doubles.sam STO14_doubles.sam
mv HTAdapter71_doubles.sam STO15_doubles.sam
mv HTAdapter72_doubles.sam STO17_doubles.sam
mv HTAdapter73_doubles.sam STO18_doubles.sam
mv HTAdapter74_doubles.sam STO19_doubles.sam
mv HTAdapter75_doubles.sam STO20_doubles.sam
mv HTAdapter76_doubles.sam STO21_doubles.sam
mv HTAdapter77_doubles.sam SUM01_doubles.sam
mv HTAdapter78_doubles.sam SUM02_doubles.sam
mv HTAdapter79_doubles.sam SUM03_doubles.sam
mv HTAdapter80_doubles.sam SUM04_doubles.sam
mv HTAdapter81_doubles.sam SUM05_doubles.sam
mv HTAdapter82_doubles.sam SUM06_doubles.sam
mv HTAdapter83_doubles.sam SUM07_doubles.sam
mv HTAdapter84_doubles.sam SUM08_doubles.sam
mv HTAdapter85_doubles.sam SUM09_doubles.sam
mv HTAdapter86_doubles.sam SUM10_doubles.sam
mv HTAdapter87_doubles.sam SUM11_doubles.sam
mv HTAdapter88_doubles.sam SUM12_doubles.sam
mv HTAdapter89_doubles.sam SUM13_doubles.sam
mv HTAdapter90_doubles.sam SUM15_doubles.sam
mv HTAdapter91_doubles.sam SUM16_doubles.sam
mv HTAdapter92_doubles.sam SUM17_doubles.sam
mv HTAdapter93_doubles.sam SUM18_doubles.sam
mv HTAdapter94_doubles.sam SUM19_doubles.sam
mv HTAdapter95_doubles.sam SUM20_doubles.sam
mv HTAdapter96_doubles.sam SUM21_doubles.sam

# create list with each sam file
 ls | grep "_doubles.sam" |cut -d "_" -f 1 | sort | uniq > readGroupList
 # while loop through the names in readGroupList
 while read i
 do
 # add read groups
 java -Xmx8g -jar /tools/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups I="$i"_doubles.sam O="$i"_IDed.sam RGID="$i"_HTAdapter          RGLB=SeqCap2012 RGPL=illumina RGPU="$i" RGSM="$i"
 done<readGroupList
 echo "RLK_report: readGrouping complete"

# make .sam files list for Samtools processing
ls | grep "_IDed.sam" |cut -d "_" -f 1 | sort | uniq  > samList
# loop through .sam files
while read i;
do
# convert .sam to .bam & sort. ex: HTAdapter84_CACCTTAC_L003_R1.sam –> HTAdapter84_CACCTTAC_L003_R1_sorted.bam
samtools view -@ 2 -bS "$i"_IDed.sam | samtools sort -@ 2 -o "$i"_sorted.bam
# index the sorted .bam
samtools index "$i"_sorted.bam
# tally mapped reads & calcuate the stats
samtools idxstats "$i"_sorted.bam > "$i"_counts.txt
samtools flagstat "$i"_sorted.bam > "$i"_stats.txt
samtools depth "$i"_sorted.bam > "$i"_depth.txt
done<samList
ls | grep "_depth.txt" | cut -d "_" -f 1 | sort | uniq > depthList
# Add individual name to each line in depth file
while read i
do
for f in "$i"_depth.txt
do
sed -i "s/$/\t$i/" $f; done
done<depthList

# Generate file with all depth information
cat *_depth.txt > Genome_depth.txt
for f in Genome_depth.txt
do
sed -i "s/$/\tGenome/" $f; done

# Create file with avg depth per exon using Randy's python script
cp /home/rlk0015/SeqCap/pythonScripts/Genome_avgDepth.py .
python Genome_avgDepth.py Genome_depth.txt Genome_avgDepth.txt

# make results directory & move results
mkdir -p /home/rlk0015/SeqCap/Dec2017Results/PE/Genome_WorkingDirectory/stats
mkdir -p /home/rlk0015/SeqCap/Dec2017Results/PE/Genome_WorkingDirectory/counts
mkdir -p /home/rlk0015/SeqCap/Dec2017Results/PE/Genome_WorkingDirectory/avgDepth

# cp results to respective directories
cp *stats.txt /home/rlk0015/SeqCap/Dec2017Results/PE/Genome_WorkingDirectory/stats
cp *counts.txt /home/rlk0015/SeqCap/Dec2017Results/PE/Genome_WorkingDirectory/counts
cp *depth.txt /home/rlk0015/SeqCap/Dec2017Results/PE/Genome_WorkingDirectory/avgDepth

