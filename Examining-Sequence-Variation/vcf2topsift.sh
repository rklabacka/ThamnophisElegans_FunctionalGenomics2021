#!/bin/sh

function filter_top_fst_dxy_vcf {
  cd $WorkingDirectory/variantFiltration
  # you should have copied top_fst_missense.bed into references
  vcftools \
    --gzvcf Full_CDS_missense.vcf.gz \
    --bed $WorkingDirectory/References/top_missense.bed \
    --out top_fst_dxy --recode --keep-INFO-all
  mv top_fst_dxy.recode.vcf top_fst_dxy.vcf
  bgzip top_fst_dxy.vcf
  bcftools index -f top_fst_dxy.vcf.gz
}

function add_sift_scores {
  cd $WorkingDirectory/variantFiltration
  bcftools isec -p sift_isec_results -n=2 Full_Exons_SIFTpredictions.vcf.gz top_fst_dxy.vcf.gz -o top_fst_dxy_sift.vcf.gz -O z
  cp sift_isec_results/0000.vcf.gz top_fst_dxy_sift.vcf.gz
}


function tab_sift_vcf {
  cd $WorkingDirectory/variantFiltration
  gunzip top_fst_dxy_sift.vcf.gz
  gunzip top_fst_dxy.vcf.gz
  # Create top_fst_dxy_sift.txt table and top_fst_dxy_sift.txt.gene_list
  python $pythonScripts/parse_sift_vcf.py top_fst_dxy_sift.vcf top_fst_dxy.vcf top_fst_dxy_sift.txt
}

function get_sift_faa {
  cd $WorkingDirectory/variantFiltration
  awk '{print $5}' top_fst_dxy_sift.txt > top_fst_dxy_sift_genes.txt
  cd $WorkingDirectory/SNP_analysis/vcf2fasta/RefSeq/Sequences/CapturedGenes
  while read i
  do
  cat *"$i"*.faa >> top_fst_dxy_sift_genes.faa
  echo >> top_fst_dxy_sift_genes.faa
  done<$WorkingDirectory/variantFiltration/top_fst_dxy_sift_genes.txt  
}

