#!/usr/bin/env Rscript

################################################################################
###Script was created by S. Sefick in 2016 to plot vcf stats in order to determine hard filtering values.
###worked with R/3.4.3


require(argparse)

parser <- ArgumentParser()
                                        # specify our desired options
                                        # by default ArgumentParser will add an help option


parser$add_argument("-I", "--input_vcf", action="store", default="Table from VCF file", dest="input_vcf",
                    help="Input %(default)s")

parser$add_argument("-O", "--output_pdf", action="store", default="out.pdf", dest="output_pdf",
                    help="Output %(default)s")

parser$add_argument("-DP", "--Depth_cut_off", action="store", default="TRUE", dest="raw_DP",
                    help="TRUE/FALSE; Truncate graph to 5*SD; defaults to  %(default)s")


                                        # get command line options, if help option encountered print help and exit,
                                        # otherwise if options not found on command line then set defaults,
argv <- parser$parse_args()

################################################################################

                                        #test

test <- 0

if(test==1){

    input <- "/home/ssefick/scratch/arctoides/ORINGINAL_GVCF_MACAQUE_COMBINED/split/chrM_Annotations/filtered_chrM_RGQ_VQSR_ANNOTATED"

    output <- paste(input, "pdf", sep=".")

    argv <- list(input_vcf=input, output_pdf=output, raw_DP=TRUE)
}

if(argv["input_vcf"][[1]]=="Table from VCF file"){stop("\n\n\n Needs 2 agruments: \n\n Argument 1: vcf table generated with GATK \n Argument2: pdf file for graphical output \n\n For full help:\n plot_stats_vcf_CL_ARGS.R --help \n\n\n")}

input <- argv["input_vcf"][[1]]
output <- argv["output_pdf"][[1]]
raw_DP <- argv["raw_DP"][[1]]

##libraries
require(ggplot2)
require(grid)
require(gtable)
require(gridExtra)
require(tidyverse)
require(dplyr)


##get data in
x <- read.table(input, sep="\t", header=TRUE, na.strings="None")
x <- x %>% filter(GQ < 20)


##A positive value would suggest ref near end; in practice only neg. values for filtering variants
##ReadPosRankSum
RPRS.plot <- ggplot(x, aes(x=ReadPosRankSum))+geom_density(alpha=0.2)
RPRS.plot <- RPRS.plot + geom_vline(xintercept = -8, colour="red")
RPRS.plot <- RPRS.plot + ggtitle("ReadPosRankSum - Negative values indicate ALT near ends")

##A positive value would suggest ref less often supporting a site; in practice only neg. values for filtering variants
##MQRankSum
MQRS.plot <- ggplot(x, aes(x=MQRankSum))+geom_density(alpha=0.2)
MQRS.plot <- MQRS.plot + geom_vline(xintercept = -12.5, colour="red")
MQRS.plot <- MQRS.plot + ggtitle("MQRankSum of Mapping Quality - Negative values indicate ALT near ends (-12.5 cutoff)")

##SOR
SOR.plot <- ggplot(x, aes(x=SOR))+geom_density(alpha=0.2)
SOR.plot <- SOR.plot + geom_vline(xintercept = 3, colour="red")
SOR.plot <- SOR.plot + ggtitle("StrandOddsRatio (SOR) - Greater than 3 shows strand bias")

##MQ
MQ.plot <- ggplot(x, aes(x=MQ))+geom_density(alpha=0.2)
MQ.plot <- MQ.plot + geom_vline(xintercept = 40, colour="red")
MQ.plot <- MQ.plot + ggtitle("Root Mean Square Error Mapping Quality (MQ) - Less than 40 removed")

##FS can be 0 start here
##plot FS data
FS.plot <- ggplot(x, aes(x=FS))+geom_density(alpha=0.2)
FS.plot <- FS.plot + geom_vline(xintercept = 60, colour="red")

FS.plot <- FS.plot + ggtitle("Fisher Strand (FS- Strand Bias Variant on F or R) - Greater than 60 removed")

##plot QD
QD.plot <- ggplot(x, aes(x=QD))+geom_density(alpha=0.2)
QD.plot <- QD.plot + geom_vline(xintercept = 2, colour="red")
QD.plot <- QD.plot + ggtitle("QD -variant qual/unfiltered depth - Below 2 removed")

##GQ
GQ.plot <- ggplot(x, aes(x=GQ))
GQ.plot <- GQ.plot + geom_density(alpha=0.2)
GQ.plot <- GQ.plot + geom_vline(xintercept = 30, colour="red")
GQ.plot <- GQ.plot + ggtitle("GQ- Genotype Quality (Phred score) - No current cutoff")

##plot DP
DP.plot <- ggplot(x, aes(x=DP))
DP.plot <- DP.plot + geom_density(alpha=0.2)
DP.plot <- DP.plot + ggtitle("DP- Depth - Red line drawn at Mean + 5*STDEV")

x_bar<- mean(x$DP, na.rm=TRUE)
SD <- sd(x$DP, na.rm=TRUE)
x_bar_SD <- x_bar + (5*SD)

stats <- data.frame(mean=x_bar, SD=SD, mean_plus_5xSD=x_bar_SD)

if(raw_DP=="TRUE"){

    DP.plot <- DP.plot + geom_vline(xintercept=x_bar_SD, col="red")+xlim(0, x_bar_SD+50)

}

pdf(output)

# FS.plot

QD.plot

RPRS.plot

MQRS.plot

SOR.plot

MQ.plot

GQ.plot

DP.plot

write.csv(x, 'snps.csv')

############################################################################
##table Depth Stats
############################################################################

##old table
##grid.arrange(tableGrob(stats))

t1 <- tableGrob(stats)
title <- textGrob("Depth Statistics",gp=gpar(fontsize=50))
padding <- unit(5,"mm")
table <- gtable_add_rows(t1, heights = grobHeight(title) + padding, pos = 0)
table <- gtable_add_grob(table, title, 1, 1, 1, ncol(table))
grid.newpage()
grid.draw(table)
############################################################################
############################################################################
############################################################################

dev.off()
