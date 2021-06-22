# MAIN
loadModules
WorkingDirectory=/scratch/rlk0015/Telag/May2020/WorkingDirectory
pythonScripts=/home/rlk0015/SeqCap/code/GenomicProcessingPipeline/pythonScripts
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

#+ COMPLETED vcf2faa
reference2faa
#+ COMPLETED moveCapturedGenes
#+ COMPLETED createMSA faa protein Sequences maskedMSA
#+ COMPLETED createMSA fna transcript Sequences maskedMSA

#+ COMPLETED createPopFiles
#+ COMPLETED createPairwiseVCFs
#+ COMPLETED getPairwisePopGen


