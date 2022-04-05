#!/bin/sh

# This is the bash script used for functional SNP annotation
# using the software SIFT

# Prepared by Randy Klabacka

#Give job a name
#PBS -N functionalAnnotation

#-- We recommend passing your environment variables down to the
#-- compute nodes with -V, but this is optional
#PBS -V

#-- Specify the number of nodes and cores you want to use
#-- Hopper's standard compute nodes have a total of 20 cores each
#-- so, to use all the processors on a single machine, set your
#-- ppn (processors per node) to 20.
#PBS -l nodes=1:ppn=20,walltime=20:00:00:00
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

WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory/SNP_analysis/proteinStructure/Thamnophis_elegans

function checkSIFT4G {
# -- Optional pre-check: see if sift4g is working
## --- I had to install sift4g locally using binutils 2.31
cd /home/rlk0015/scripts_to_build_SIFT_db
perl make-SIFT-db-all.pl -config test_files/candidatus_carsonella_ruddii_pv_config.txt --ensembl_download
}

function createWorkingEnvironment-sift2vcf {
# -- Step 1: Set up working environment based on https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB
mkdir -p $WorkingDirectory
  cd $WorkingDirectory
    mkdir -p gene-annotation-src
    mkdir -p chr-src
    mkdir -p dbSNP
    mkdir -p proteinDB
}

function copyRef-sift2vcf {
# -- Step 2: Copy reference genome and annotation into working environment
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/References
cp TelagGenome.fasta Thamnophis_elegans.fa
cp TelagGenome.gtf Thamnophis_elegans.gtf
bgzip Thamnophis_elegans.fa
bgzip Thamnophis_elegans.gtf
mv Thamnophis_elegans.fa.gz $WorkingDirectory/chr-src
mv Thamnophis_elegans.gtf.gz $WorkingDirectory/gene-annotation-src
}

function downloadUniRef {
# -- Step 3: Download protein database (this will take a long time)
cd $WorkingDirectory/proteinDB
wget ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
}

# -- Step 4: Create configuration file (Thamnophis_elegans.txt) following directions in
# -- https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB

function createSIFTdatabase {
# -- Step 5: Create database
## --- I had an issue with my alignment while creating my database.
## --- I followed rvaser's advice here: https://github.com/rvaser/sift4g/issues/10
## --- This should only be done if issues with alignment (specifically the error 
## --- defined by jarobin) are happening.
cd /home/rlk0015/scripts_to_build_SIFT_db
/tools/perl-5.26.1/bin/perl make-SIFT-db-all.pl -config $WorkingDirectory/Thamnophis_elegans.txt
}

function annotateVCF-sift2vcf {
# -- Step 6: Annotate vcf with SIFT scores
## --- This is done with the SIFT4G_Annotator
## --- Clone from here: https://github.com/pauline-ng/SIFT4G_Annotator
cd /scratch/rlk0015/Telag/May2020/WorkingDirectory/variantFiltration
cp Full_Exons.vcf.gz $WorkingDirectory/dbSNP/
cd $WorkingDirectory/dbSNP
gunzip Full_Exons.vcf.gz
cd ..
mkdir -p SIFT_annotation
java -jar /home/rlk0015/sift4g_annotator/SIFT4G_Annotator.jar -c -i $WorkingDirectory/dbSNP/Full_Exons.vcf -d $WorkingDirectory/rThaEle1.pri -r $WorkingDirectory/SIFT_annotation -t
}
