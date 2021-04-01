library(tidyverse)
library(ggplot2)
library(RColorBrewer)

setwd("~/Desktop/garter_snakes/cleaned up/")

# Per Gene Meadow-v-Lake Shell Output -----------------------------------------------------------
fst=read.csv("per_gene_shell_outputs/genes.meadow-v-lake.csv",header=T,na.strings = "NaN")
categories=read.csv("SeqCapGeneCategories.csv",header=T)
fst$gene2=str_remove(fst$gene,"gene-")
categories$Category_I=str_replace(categories$Category_I,"oxidative_phosporylation","oxphos")
categories$Category_I=str_replace(categories$Category_I,"metabolsim","metabolism")
categories$Category_I=str_replace(categories$Category_I,"Stress","stress")

merged=merge(fst,categories,by.x="gene2",by.y="GarterSnakeGenome_GeneID",all=T)
write.csv(merged,file="per_gene_R_outputs/pop_gen_summary_lake_v_meadow.csv",quote=F,row.names = F)

summary(merged$Fst_Meadow_Lake)
hist(merged$Fst_Meadow_Lake)
d=density(merged$Fst_Meadow_Lake,na.rm=T)
plot(d,main="Fst density")
polygon(d,col="blue")
#abline(v=0.04,col="red")

summary(merged$dxy_Meadow_Lake)
hist(merged$dxy_Meadow_Lake)
d=density(merged$dxy_Meadow_Lake,na.rm=T)
plot(d,main="Dxy density")
polygon(d,col="blue")
#abline(v=0.345,col="red")

# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

merged$Category_I=as.factor(merged$Category_I)

pdf(file="per_gene_R_outputs/pop_gen_stats_summary.pdf")

fst.plot <- ggplot(merged, aes(x=Fst_Meadow_Lake,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
fst.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")


dxy.plot <- ggplot(merged, aes(x=dxy_Meadow_Lake,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
dxy.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")


pi1.plot <- ggplot(merged, aes(x=pi_Meadow,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
pi1.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")


pi2.plot <- ggplot(merged, aes(x=pi_Lake,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
pi2.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")

dev.off()

fst.mean=as.data.frame(tapply(merged$Fst_Meadow_Lake,merged$Category_I,mean,na.rm=T))
dxy.mean=as.data.frame(tapply(merged$dxy_Meadow_Lake,merged$Category_I,mean,na.rm=T))


#plot fst vs. dxy
plot(merged$Fst_Meadow_Lake,merged$dxy_Meadow_Lake)
cor.test(merged$Fst_Meadow_Lake,merged$dxy_Meadow_Lake)

#fst.top=subset(merged,merged$Fst_Meadow_Lake>quantile(merged$Fst_Meadow_Lake,na.rm=T)[4])
#dxy.top=subset(merged,merged$dxy_Meadow_Lake>quantile(merged$dxy_Meadow_Lake,na.rm=T)[4])

#top.merged=merge(fst.top,dxy.top,by="gene2")
#write.table(top.merged[c(1,9:10)],file="top.dxy.fst.tsv",quote=F,sep="\t",row.names = F)



# Per Gene All Populations ---------------------------------------------------------

all_fst=read.csv("per_gene_shell_outputs/genes.all-pops.csv",header=T,na.strings = "NaN")
#categories=read.csv("SeqCapGeneCategories.csv",header=T)
all_fst$gene2=str_remove(all_fst$gene,"gene-")
#categories$Category_I=str_replace(categories$Category_I,"oxidative_phosporylation","oxphos")
#categories$Category_I=str_replace(categories$Category_I,"metabolsim","metabolism")
#categories$Category_I=str_replace(categories$Category_I,"Stress","stress")

new_merged=merge(all_fst,categories,by.x="gene2",by.y="GarterSnakeGenome_GeneID",all=T)
write.csv(new_merged,file="per_gene_R_outputs/pop_gen_summary_all-pops.csv",quote=F,row.names = F)



# Per Gene Fst for pairwise lake vs meadow comparison to lake-v-meadow estimates --------

all_merged=merge(all_fst,fst,by="gene2",all=T)

#all_merged=subset(all_merged,!is.nan(all_merged$Fst_Meadow_Lake))

#MAH Meadow
#MER Lake
#PVM Meadow
#SUM Meadow
#STO Lake
#CHR Lake
#RON Meadow
#ROC Lake
#ELF Lake
#NAM Meadow
#MAR Lake
#PAP Meadow

within=vector(mode="numeric",length=length(all_merged$gene2))
between=vector(mode="numeric",length=length(all_merged$gene2))
mean_l_v_l=vector(mode="numeric",length=length(all_merged$gene2))
mean_l_v_m=vector(mode="numeric",length=length(all_merged$gene2))
mean_m_v_m=vector(mode="numeric",length=length(all_merged$gene2))
mean_diff=vector(mode="numeric",length=length(all_merged$gene2))


for (g in 1:length(all_merged$gene2)) {
  l_v_l=c(all_merged$Fst_MER_STO[g],all_merged$Fst_MER_CHR[g],all_merged$Fst_MER_ROC[g],
          all_merged$Fst_MER_ELF[g],all_merged$Fst_MER_MAR[g],all_merged$Fst_STO_CHR[g],
          all_merged$Fst_STO_ROC[g],all_merged$Fst_STO_ELF[g],all_merged$Fst_STO_MAR[g],
          all_merged$Fst_CHR_ROC[g],all_merged$Fst_CHR_ELF[g],all_merged$Fst_CHR_MAR[g],
          all_merged$Fst_ROC_ELF[g],all_merged$Fst_ROC_MAR[g],all_merged$Fst_ELF_MAR[g])
  m_v_m=c(all_merged$Fst_MAH_PVM[g],all_merged$Fst_MAH_SUM[g],all_merged$Fst_MAH_RON[g],
          all_merged$Fst_MAH_NAM[g],all_merged$Fst_MAH_PAP[g],all_merged$Fst_PVM_SUM[g],
          all_merged$Fst_PVM_RON[g],all_merged$Fst_PVM_NAM[g],all_merged$Fst_PVM_PAP[g],
          all_merged$Fst_SUM_RON[g],all_merged$Fst_SUM_NAM[g],all_merged$Fst_SUM_PAP[g],
          all_merged$Fst_RON_NAM[g],all_merged$Fst_RON_PAP[g],all_merged$Fst_NAM_PAP[g])
  l_v_m=c(all_merged$Fst_MAH_MER[g],all_merged$Fst_MAH_STO[g],all_merged$Fst_MAH_CHR[g],
          all_merged$Fst_MAH_ROC[g],all_merged$Fst_MAH_ELF[g],all_merged$Fst_MAH_MAR[g],
          all_merged$Fst_MER_PVM[g],all_merged$Fst_MER_SUM[g],all_merged$Fst_MER_RON[g],
          all_merged$Fst_MER_NAM[g],all_merged$Fst_MER_PAP[g],all_merged$Fst_PVM_STO[g],
          all_merged$Fst_PVM_CHR[g],all_merged$Fst_PVM_ROC[g],all_merged$Fst_PVM_ELF[g],
          all_merged$Fst_PVM_MAR[g],all_merged$Fst_SUM_STO[g],all_merged$Fst_SUM_CHR[g],
          all_merged$Fst_SUM_ROC[g],all_merged$Fst_SUM_ELF[g],all_merged$Fst_SUM_MAR[g],
          all_merged$Fst_STO_RON[g],all_merged$Fst_STO_NAM[g],all_merged$Fst_STO_PAP[g],
          all_merged$Fst_CHR_RON[g],all_merged$Fst_CHR_NAM[g],all_merged$Fst_CHR_PAP[g],
          all_merged$Fst_RON_ROC[g],all_merged$Fst_RON_ELF[g],all_merged$Fst_RON_MAR[g],
          all_merged$Fst_ROC_NAM[g],all_merged$Fst_ROC_PAP[g],all_merged$Fst_ELF_NAM[g],
          all_merged$Fst_ELF_PAP[g],all_merged$Fst_NAM_MAR[g],all_merged$Fst_MAR_PAP[g])
  
  mean_l_v_l[g]=mean(l_v_l,na.rm=T)
  mean_m_v_m[g]=mean(m_v_m,na.rm=T)
  mean_l_v_m[g]=mean(l_v_m,na.rm=T)
  
  #combine l_v_l and m_v_m to do ttest below
  combined=c(l_v_l,m_v_m)
  
  #skip 
  if (all(is.na(l_v_l)) || all(is.na(m_v_m)) || all(is.na(l_v_m)) || sum(!is.na(l_v_l))==1 || sum(!is.na(l_v_m))==1 || sum(!is.na(m_v_m))==1) { 
    skip=all_merged$gene[g]
    print(skip)
    between[g]=NaN
    within[g]=NaN
    mean_diff[g]=NaN
  } else {
    x=t.test(l_v_l,m_v_m)
    within[g]=x$p.value
    y=t.test(combined,l_v_m)
    between[g]=y$p.value
    mean_diff[g]=mean(l_v_m,na.rm=T)-mean(combined,na.rm=T)
    #fit <- lm(MisSynProp~EcotypeComp, data = dat)
    #slope <- coef(fit)[2]
    #pVal <- anova(fit)$'Pr(>F)'[1]
    
  }
  
}

all_merged$within=within
all_merged$between=between
all_merged$l_v_l_mean=mean_l_v_l
all_merged$l_v_m_mean=mean_l_v_m
all_merged$m_v_m_mean=mean_m_v_m
all_merged$mean_diff=mean_diff

fst_gene_list=all_merged$gene2[all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean & all_merged$mean_diff>0]
unique(fst_gene_list[!is.na(fst_gene_list)])

#top_fst=subset(all_merged,all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean)
top_fst=all_merged[,c(1:4,160,162:167)]

write.csv(top_fst,file="per_gene_R_outputs/fst_within_v_between.csv",quote=F,row.names = F)

fst_final=merge(top_fst,categories,by.x="gene2",by.y="GarterSnakeGenome_GeneID",all=T)
signif_fst=subset(fst_final,fst_final$between<0.05 & fst_final$l_v_m_mean>fst_final$l_v_l_mean & fst_final$l_v_m_mean>fst_final$m_v_m_mean & fst_final$mean_diff>0)



# Per Gene Density Plot of FST results ---------------------------------------------

#combined with categories
fst_final$Category_I=as.factor(fst_final$Category_I)
plot(fst_final$mean_diff~fst_final$Category_I,main="Mean difference in average Fst within ecotypes versus between ecotypes",xlab="Gene Categories",ylab="Mean difference in FST values", type="n")
abline(h=0, col="grey",lty=2 )
plot(fst_final$mean_diff~fst_final$Category_I,add=T, col=rainbow(14))
boxplot(signif_fst$mean_diff~signif_fst$Category_I,col=rainbow(14),add=T)

#look at fst distribution in these genes
behavior=density(fst_final$mean_diff[fst_final$Category_I=="behavior"],na.rm = T)
color=density(fst_final$mean_diff[fst_final$Category_I=="color"],na.rm = T)
growth=density(fst_final$mean_diff[fst_final$Category_I=="growth"],na.rm = T)
hypoxia=density(fst_final$mean_diff[fst_final$Category_I=="hypoxia"],na.rm = T)
IILS=density(fst_final$mean_diff[fst_final$Category_I=="IILS"],na.rm = T)
immune=density(fst_final$mean_diff[fst_final$Category_I=="immune"],na.rm = T)
metabolism=density(fst_final$mean_diff[fst_final$Category_I=="metabolism"],na.rm = T)
oxidative_stress=density(fst_final$mean_diff[fst_final$Category_I=="oxidative_stress"],na.rm = T)
oxphos=density(fst_final$mean_diff[fst_final$Category_I=="oxphos"],na.rm = T)
p53=density(fst_final$mean_diff[fst_final$Category_I=="p53"],na.rm = T)
PI3K=density(fst_final$mean_diff[fst_final$Category_I=="PI3K"],na.rm = T)
random=density(fst_final$mean_diff[fst_final$Category_I=="random"],na.rm = T)
reproduction=density(fst_final$mean_diff[fst_final$Category_I=="reproduction"],na.rm = T)
stress=density(fst_final$mean_diff[fst_final$Category_I=="stress"],na.rm = T)

#get max x/y values
x_max=max(c(behavior$x,color$x,growth$x,hypoxia$x,IILS$x,immune$x,metabolism$x,oxidative_stress$x,oxphos$x,p53$x,PI3K$x,random$x,reproduction$x,stress$x))
x_min=min(c(behavior$x,color$x,growth$x,hypoxia$x,IILS$x,immune$x,metabolism$x,oxidative_stress$x,oxphos$x,p53$x,PI3K$x,random$x,reproduction$x,stress$x))
y_max=max(c(behavior$y,color$y,growth$y,hypoxia$y,IILS$y,immune$y,metabolism$y,oxidative_stress$y,oxphos$y,p53$y,PI3K$y,random$y,reproduction$y,stress$y))
y_min=min(c(behavior$y,color$y,growth$y,hypoxia$y,IILS$y,immune$y,metabolism$y,oxidative_stress$y,oxphos$y,p53$y,PI3K$y,random$y,reproduction$y,stress$y))

pdf(file="per_gene_R_outputs/mean_diff_fst.pdf",width=8.5,height=11)
mycols=colorRampPalette(brewer.pal(9,"Set1"))(14)
#d=density(top_fst$Fst_Meadow_Lake,na.rm=T)
plot(behavior,xlim=c(x_min,x_max),ylim=c(y_min,y_max),type="n",main="Density plots of mean difference of FST values within vs between ecotypes",xlab="Mean difference in FST values")
lines(behavior,col=mycols[1],lwd=2)
lines(color,col=mycols[2],lwd=2)
lines(growth,col=mycols[3],lwd=2)
lines(hypoxia,col=mycols[4],lwd=2)
lines(IILS,col=mycols[5],lwd=2)
lines(immune,col=mycols[6],lwd=2)
lines(metabolism,col=mycols[7],lwd=2)
lines(oxidative_stress,col=mycols[8],lwd=2)
lines(oxphos,col=mycols[9],lwd=2)
lines(p53,col=mycols[10],lwd=2)
lines(PI3K,col=mycols[11],lwd=2)
lines(random,col=mycols[12],lwd=2)
lines(reproduction,col=mycols[13],lwd=2)
lines(stress,col=mycols[14],lwd=2)
#abline(v=0.04,col="red")
legend("topright",legend=c("behavior","color","growth","hypoxia","IILS","immune",
                           "metabolism","oxidative_stress","oxphos","p53",
                           "PI3K","random","reproduction","stress"),col=mycols,
       lty=1,lwd=2,cex=1.5,bty="n")
#text(locator(), labels = c("color","IILS","metabolism"))
dev.off()

#make same plot but ONLY with significant genes
#behavior=density(signif_fst$mean_diff[signif_fst$Category_I=="behavior"],na.rm = T)
#color=density(signif_fst$mean_diff[signif_fst$Category_I=="color"],na.rm = T)
#growth=density(signif_fst$mean_diff[signif_fst$Category_I=="growth"],na.rm = T)
#hypoxia=density(signif_fst$mean_diff[signif_fst$Category_I=="hypoxia"],na.rm = T)
#IILS=density(signif_fst$mean_diff[signif_fst$Category_I=="IILS"],na.rm = T)
#immune=density(signif_fst$mean_diff[signif_fst$Category_I=="immune"],na.rm = T)
#metabolism=density(signif_fst$mean_diff[signif_fst$Category_I=="metabolism"],na.rm = T)
#oxidative_stress=density(signif_fst$mean_diff[signif_fst$Category_I=="oxidative_stress"],na.rm = T)
#oxphos=density(signif_fst$mean_diff[signif_fst$Category_I=="oxphos"],na.rm = T)
#p53=density(signif_fst$mean_diff[signif_fst$Category_I=="p53"],na.rm = T)
#PI3K=density(signif_fst$mean_diff[signif_fst$Category_I=="PI3K"],na.rm = T)
#random=density(signif_fst$mean_diff[signif_fst$Category_I=="random"],na.rm = T)
#reproduction=density(signif_fst$mean_diff[signif_fst$Category_I=="reproduction"],na.rm = T)
#stress=density(signif_fst$mean_diff[signif_fst$Category_I=="stress"],na.rm = T)

#get max x/y values
#x_max=max(c(behavior$x,color$x,growth$x,hypoxia$x,IILS$x,immune$x,metabolism$x,oxidative_stress$x,oxphos$x,p53$x,PI3K$x,random$x,reproduction$x,stress$x))
#x_min=min(c(behavior$x,color$x,growth$x,hypoxia$x,IILS$x,immune$x,metabolism$x,oxidative_stress$x,oxphos$x,p53$x,PI3K$x,random$x,reproduction$x,stress$x))
#y_max=max(c(behavior$y,color$y,growth$y,hypoxia$y,IILS$y,immune$y,metabolism$y,oxidative_stress$y,oxphos$y,p53$y,PI3K$y,random$y,reproduction$y,stress$y))
#y_min=min(c(behavior$y,color$y,growth$y,hypoxia$y,IILS$y,immune$y,metabolism$y,oxidative_stress$y,oxphos$y,p53$y,PI3K$y,random$y,reproduction$y,stress$y))

#pdf(file="per_gene_R_outputs/mean_diff_fst_signif_ONLY.pdf",width=8.5,height=11)
#mycols=colorRampPalette(brewer.pal(9,"Set1"))(14)
#d=density(top_fst$Fst_Meadow_Lake,na.rm=T)
#plot(behavior,xlim=c(x_min,x_max),ylim=c(y_min,y_max),type="n",main="Density plots of mean difference of FST values within vs between ecotypes",xlab="Mean difference in FST values")
#lines(behavior,col=mycols[1],lwd=2)
#lines(color,col=mycols[2],lwd=2)
#lines(growth,col=mycols[3],lwd=2)
#lines(hypoxia,col=mycols[4],lwd=2)
#lines(IILS,col=mycols[5],lwd=2)
#lines(immune,col=mycols[6],lwd=2)
#lines(metabolism,col=mycols[7],lwd=2)
#lines(oxidative_stress,col=mycols[8],lwd=2)
#lines(oxphos,col=mycols[9],lwd=2)
#lines(p53,col=mycols[10],lwd=2)
#lines(PI3K,col=mycols[11],lwd=2)
#lines(random,col=mycols[12],lwd=2)
#lines(reproduction,col=mycols[13],lwd=2)
#lines(stress,col=mycols[14],lwd=2)
#abline(v=0.04,col="red")
#legend("topright",legend=c("behavior","color","growth","hypoxia","IILS","immune",
#                           "metabolism","oxidative_stress","oxphos","p53",
#                           "PI3K","random","reproduction","stress"),col=mycols,
#       lty=1,lwd=2,cex=1.5,bty="n")
##text(locator(), labels = c("color","IILS","metabolism"))
#dev.off()


# Per Gene Repeat, but with DXY ----------------------------------------------------

#repeat, but with DXY
all_merged=merge(all_fst,fst,by="gene2",all=T)

#all_merged=subset(all_merged,!is.nan(all_merged$dxy_Meadow_Lake))

within=vector(mode="numeric",length=length(all_merged$gene2))
between=vector(mode="numeric",length=length(all_merged$gene2))
mean_l_v_l=vector(mode="numeric",length=length(all_merged$gene2))
mean_l_v_m=vector(mode="numeric",length=length(all_merged$gene2))
mean_m_v_m=vector(mode="numeric",length=length(all_merged$gene2))
mean_diff=vector(mode="numeric",length=length(all_merged$gene2))


for (g in 1:length(all_merged$gene2)) {
  l_v_l=c(all_merged$dxy_MER_STO[g],all_merged$dxy_MER_CHR[g],all_merged$dxy_MER_ROC[g],
          all_merged$dxy_MER_ELF[g],all_merged$dxy_MER_MAR[g],all_merged$dxy_STO_CHR[g],
          all_merged$dxy_STO_ROC[g],all_merged$dxy_STO_ELF[g],all_merged$dxy_STO_MAR[g],
          all_merged$dxy_CHR_ROC[g],all_merged$dxy_CHR_ELF[g],all_merged$dxy_CHR_MAR[g],
          all_merged$dxy_ROC_ELF[g],all_merged$dxy_ROC_MAR[g],all_merged$dxy_ELF_MAR[g])
  m_v_m=c(all_merged$dxy_MAH_PVM[g],all_merged$dxy_MAH_SUM[g],all_merged$dxy_MAH_RON[g],
          all_merged$dxy_MAH_NAM[g],all_merged$dxy_MAH_PAP[g],all_merged$dxy_PVM_SUM[g],
          all_merged$dxy_PVM_RON[g],all_merged$dxy_PVM_NAM[g],all_merged$dxy_PVM_PAP[g],
          all_merged$dxy_SUM_RON[g],all_merged$dxy_SUM_NAM[g],all_merged$dxy_SUM_PAP[g],
          all_merged$dxy_RON_NAM[g],all_merged$dxy_RON_PAP[g],all_merged$dxy_NAM_PAP[g])
  l_v_m=c(all_merged$dxy_MAH_MER[g],all_merged$dxy_MAH_STO[g],all_merged$dxy_MAH_CHR[g],
          all_merged$dxy_MAH_ROC[g],all_merged$dxy_MAH_ELF[g],all_merged$dxy_MAH_MAR[g],
          all_merged$dxy_MER_PVM[g],all_merged$dxy_MER_SUM[g],all_merged$dxy_MER_RON[g],
          all_merged$dxy_MER_NAM[g],all_merged$dxy_MER_PAP[g],all_merged$dxy_PVM_STO[g],
          all_merged$dxy_PVM_CHR[g],all_merged$dxy_PVM_ROC[g],all_merged$dxy_PVM_ELF[g],
          all_merged$dxy_PVM_MAR[g],all_merged$dxy_SUM_STO[g],all_merged$dxy_SUM_CHR[g],
          all_merged$dxy_SUM_ROC[g],all_merged$dxy_SUM_ELF[g],all_merged$dxy_SUM_MAR[g],
          all_merged$dxy_STO_RON[g],all_merged$dxy_STO_NAM[g],all_merged$dxy_STO_PAP[g],
          all_merged$dxy_CHR_RON[g],all_merged$dxy_CHR_NAM[g],all_merged$dxy_CHR_PAP[g],
          all_merged$dxy_RON_ROC[g],all_merged$dxy_RON_ELF[g],all_merged$dxy_RON_MAR[g],
          all_merged$dxy_ROC_NAM[g],all_merged$dxy_ROC_PAP[g],all_merged$dxy_ELF_NAM[g],
          all_merged$dxy_ELF_PAP[g],all_merged$dxy_NAM_MAR[g],all_merged$dxy_MAR_PAP[g])
  
  mean_l_v_l[g]=mean(l_v_l,na.rm=T)
  mean_m_v_m[g]=mean(m_v_m,na.rm=T)
  mean_l_v_m[g]=mean(l_v_m,na.rm=T)
  
  #combine l_v_l and m_v_m to do ttest below
  combined=c(l_v_l,m_v_m)
  
  #skip 
  if (all(is.na(l_v_l)) || all(is.na(m_v_m)) || all(is.na(l_v_m)) || sum(!is.na(l_v_l))==1 || sum(!is.na(l_v_m))==1 || sum(!is.na(m_v_m))==1) { 
    skip=all_merged$gene[g]
    print(skip)
    between[g]=NaN
    within[g]=NaN
    mean_diff[g]=NaN
  } else {
    x=t.test(l_v_l,m_v_m)
    within[g]=x$p.value
    y=t.test(combined,l_v_m)
    between[g]=y$p.value
    mean_diff[g]=mean(l_v_m,na.rm=T)-mean(combined,na.rm=T)
    #fit <- lm(MisSynProp~EcotypeComp, data = dat)
    #slope <- coef(fit)[2]
    #pVal <- anova(fit)$'Pr(>F)'[1]
    
  }
  
}

all_merged$within=within
all_merged$between=between
all_merged$l_v_l_mean=mean_l_v_l
all_merged$l_v_m_mean=mean_l_v_m
all_merged$m_v_m_mean=mean_m_v_m
all_merged$mean_diff=mean_diff

dxy_gene_list=all_merged$gene2[all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean & all_merged$mean_diff>0]
dxy_gene_list[!is.na(dxy_gene_list)]

#top_dxy=subset(all_merged,all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean)
top_dxy=all_merged[,c(1:4,159,162:167)]
write.csv(top_dxy,file="per_gene_R_outputs/dxy_within_v_between.csv",quote=F,row.names = F)


#intersect fst/dxy lists
intersect(fst_gene_list[!is.na(fst_gene_list)],dxy_gene_list[!is.na(dxy_gene_list)])

dxy_final=merge(top_dxy,categories,by.x="gene2",by.y="GarterSnakeGenome_GeneID",all=T)
signif_dxy=subset(dxy_final,dxy_final$between<0.05 & dxy_final$l_v_m_mean>dxy_final$l_v_l_mean & dxy_final$l_v_m_mean>dxy_final$m_v_m_mean & dxy_final$mean_diff>0)



# Per Gene Density Plot of DXY results ---------------------------------------------

#combined with categories
dxy_final$Category_I=as.factor(dxy_final$Category_I)
plot(dxy_final$mean_diff~dxy_final$Category_I,main="Mean difference in average Fst within ecotypes versus between ecotypes",xlab="Gene Categories",ylab="Mean difference in FST values", type="n")
abline(h=0, col="grey",lty=2 )
plot(dxy_final$mean_diff~dxy_final$Category_I,add=T, col=rainbow(14))
boxplot(signif_dxy$mean_diff~signif_dxy$Category_I,col=rainbow(14),na.rm=T)

#look at fst distribution in these genes
behavior=density(dxy_final$mean_diff[dxy_final$Category_I=="behavior"],na.rm = T)
color=density(dxy_final$mean_diff[dxy_final$Category_I=="color"],na.rm = T)
growth=density(dxy_final$mean_diff[dxy_final$Category_I=="growth"],na.rm = T)
hypoxia=density(dxy_final$mean_diff[dxy_final$Category_I=="hypoxia"],na.rm = T)
IILS=density(dxy_final$mean_diff[dxy_final$Category_I=="IILS"],na.rm = T)
immune=density(dxy_final$mean_diff[dxy_final$Category_I=="immune"],na.rm = T)
metabolism=density(dxy_final$mean_diff[dxy_final$Category_I=="metabolism"],na.rm = T)
oxidative_stress=density(dxy_final$mean_diff[dxy_final$Category_I=="oxidative_stress"],na.rm = T)
oxphos=density(dxy_final$mean_diff[dxy_final$Category_I=="oxphos"],na.rm = T)
p53=density(dxy_final$mean_diff[dxy_final$Category_I=="p53"],na.rm = T)
PI3K=density(dxy_final$mean_diff[dxy_final$Category_I=="PI3K"],na.rm = T)
random=density(dxy_final$mean_diff[dxy_final$Category_I=="random"],na.rm = T)
reproduction=density(dxy_final$mean_diff[dxy_final$Category_I=="reproduction"],na.rm = T)
stress=density(dxy_final$mean_diff[dxy_final$Category_I=="stress"],na.rm = T)

#get max x/y values
x_max=max(c(behavior$x,color$x,growth$x,hypoxia$x,IILS$x,immune$x,metabolism$x,oxidative_stress$x,oxphos$x,p53$x,PI3K$x,random$x,reproduction$x,stress$x))
x_min=min(c(behavior$x,color$x,growth$x,hypoxia$x,IILS$x,immune$x,metabolism$x,oxidative_stress$x,oxphos$x,p53$x,PI3K$x,random$x,reproduction$x,stress$x))
y_max=max(c(behavior$y,color$y,growth$y,hypoxia$y,IILS$y,immune$y,metabolism$y,oxidative_stress$y,oxphos$y,p53$y,PI3K$y,random$y,reproduction$y,stress$y))
y_min=min(c(behavior$y,color$y,growth$y,hypoxia$y,IILS$y,immune$y,metabolism$y,oxidative_stress$y,oxphos$y,p53$y,PI3K$y,random$y,reproduction$y,stress$y))

pdf(file="per_gene_R_outputs/mean_diff_dxy.pdf",width=8.5,height=11)
mycols=colorRampPalette(brewer.pal(9,"Set1"))(14)
#d=density(top_fst$Fst_Meadow_Lake,na.rm=T)
plot(behavior,xlim=c(x_min,x_max),ylim=c(y_min,y_max),type="n",main="Density plots of mean difference of DXY values within vs between ecotypes",xlab="Mean difference in DXY values")
lines(behavior,col=mycols[1],lwd=2)
lines(color,col=mycols[2],lwd=2)
lines(growth,col=mycols[3],lwd=2)
lines(hypoxia,col=mycols[4],lwd=2)
lines(IILS,col=mycols[5],lwd=2)
lines(immune,col=mycols[6],lwd=2)
lines(metabolism,col=mycols[7],lwd=2)
lines(oxidative_stress,col=mycols[8],lwd=2)
lines(oxphos,col=mycols[9],lwd=2)
lines(p53,col=mycols[10],lwd=2)
lines(PI3K,col=mycols[11],lwd=2)
lines(random,col=mycols[12],lwd=2)
lines(reproduction,col=mycols[13],lwd=2)
lines(stress,col=mycols[14],lwd=2)
#abline(v=0.04,col="red")
legend("topright",legend=c("behavior","color","growth","hypoxia","IILS","immune",
                           "metabolism","oxidative_stress","oxphos","p53",
                           "PI3K","random","reproduction","stress"),col=mycols,
       lty=1,lwd=2,cex=1.5,bty="n")
#text(locator(), labels = c("color","IILS","metabolism"))
dev.off()

#make same plot but ONLY with significant genes
#hypoxia=density(signif_dxy$mean_diff[signif_dxy$Category_I=="hypoxia"],na.rm = T)
#p53=density(signif_dxy$mean_diff[signif_dxy$Category_I=="p53"],na.rm = T)
#PI3K=density(signif_dxy$mean_diff[signif_dxy$Category_I=="PI3K"],na.rm = T)
#random=density(signif_dxy$mean_diff[signif_dxy$Category_I=="random"],na.rm = T)
#stress=density(signif_dxy$mean_diff[signif_dxy$Category_I=="stress"],na.rm = T)

#get max x/y values
#x_max=max(c(hypoxia$x,p53$x,PI3K$x,random$x,stress$x))
#x_min=min(c(hypoxia$x,p53$x,PI3K$x,random$x,stress$x))
#y_max=max(c(hypoxia$y,p53$y,PI3K$y,random$y,stress$y))
#y_min=min(c(hypoxia$y,p53$y,PI3K$y,random$y,stress$y))

#pdf(file="per_gene_R_outputs/mean_diff_dxy_signif_ONLY.pdf",width=8.5,height=11)
#mycols=colorRampPalette(brewer.pal(9,"Set1"))(14)
##d=density(top_fst$Fst_Meadow_Lake,na.rm=T)
#plot(hypoxia,xlim=c(x_min,x_max),ylim=c(y_min,y_max),type="n",main="Density plots of mean difference of DXY values within vs between ecotypes",xlab="Mean difference in DXY values")
#lines(hypoxia,col=mycols[4],lwd=2)
#lines(p53,col=mycols[10],lwd=2)
#lines(PI3K,col=mycols[11],lwd=2)
#lines(random,col=mycols[12],lwd=2)
#lines(stress,col=mycols[14],lwd=2)
####abline(v=0.04,col="red")
#legend("topright",legend=c("hypoxia","p53","PI3K","random","stress"),col=mycols[c(4,10:12,14)],
#       lty=1,lwd=2,cex=1.5,bty="n")
####text(locator(), labels = c("color","IILS","metabolism"))
#dev.off()

# Per Site Missense Variants Only Meadow-v-Lake --------------------------------------------------

#repeat, but with ONLY missense variants

####setwd("~/Desktop/garter_snakes/cleaned up/per_site_shell_outputs/")

#read in data
fst=read.csv("per_site_shell_outputs/missense-sites.meadow-v-lake.csv",header=T)
#fst=read.csv("missense-sites.all-pops.csv",header=T)

#categories file
categories=read.csv("SeqCapGeneCategories.csv",header=T)
#fst$gene2=str_remove(fst$gene,"gene-")
categories$Category_I=str_replace(categories$Category_I,"oxidative_phosporylation","oxphos")
categories$Category_I=str_replace(categories$Category_I,"metabolsim","metabolism")
categories$Category_I=str_replace(categories$Category_I,"Stress","stress")

merged=merge(fst,categories,by.x="INFO",by.y="GarterSnakeGenome_GeneID",all=T)
write.csv(merged,file="per_site_R_outputs/pop_gen_summary_lake_v_meadow.missense-only.csv",quote=F,row.names = F)


summary(merged$Fst_Meadow_Lake)
hist(merged$Fst_Meadow_Lake)
d=density(merged$Fst_Meadow_Lake,na.rm=T)
plot(d,main="Fst density")
polygon(d,col="blue")


# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

merged$Category_I=as.factor(merged$Category_I)

pdf(file="per_site_R_outputs/pop_gen_stats_summary.missense-only.pdf")

fst.plot <- ggplot(merged, aes(x=Fst_Meadow_Lake,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
fst.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")


dxy.plot <- ggplot(merged, aes(x=dxy_Meadow_Lake,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
dxy.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")


pi1.plot <- ggplot(merged, aes(x=pi_Meadow,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
pi1.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")


pi2.plot <- ggplot(merged, aes(x=pi_Lake,y=Category_I,color=Category_I)) + geom_violin()+coord_flip()
pi2.plot + stat_summary(fun.data=data_summary, mult=1,geom="pointrange", color="red")

dev.off()

fst.mean=as.data.frame(tapply(merged$Fst_Meadow_Lake,merged$Category_I,mean,na.rm=T))
dxy.mean=as.data.frame(tapply(merged$dxy_Meadow_Lake,merged$Category_I,mean,na.rm=T))


#plot fst vs. dxy
plot(merged$Fst_Meadow_Lake,merged$dxy_Meadow_Lake)
cor.test(merged$Fst_Meadow_Lake,merged$dxy_Meadow_Lake)

#fst.top=subset(merged,merged$Fst_Meadow_Lake>quantile(merged$Fst_Meadow_Lake,na.rm=T)[4])
#dxy.top=subset(merged,merged$dxy_Meadow_Lake>quantile(merged$dxy_Meadow_Lake,na.rm=T)[4])

#top.merged=merge(fst.top,dxy.top,by="INFO")
#write.table(top.merged[c(1,9:10)],file="per_site_R_outputs/top.dxy.fst.missense-only.tsv",quote=F,sep="\t",row.names = F)



# Per Site All Pops Missense Variants ----------------------------------------------

all_fst=read.csv("per_site_shell_outputs/missense-sites.all-pops.csv",header=T)
categories=read.csv("SeqCapGeneCategories.csv",header=T)
#fst$gene2=str_remove(fst$gene,"gene-")
categories$Category_I=str_replace(categories$Category_I,"oxidative_phosporylation","oxphos")
categories$Category_I=str_replace(categories$Category_I,"metabolsim","metabolism")
categories$Category_I=str_replace(categories$Category_I,"Stress","stress")

merged=merge(all_fst,categories,by.x="INFO",by.y="GarterSnakeGenome_GeneID",all=T)
write.csv(merged,file="per_site_R_outputs/pop_gen_summary_all-pops.missense-only.csv",quote=F,row.names = F)


# Per Site Pairwise fst/dxy for missense sites -----------------------------------------

all_merged=merge(all_fst,fst,by="start",all=T)

#all_merged=subset(all_merged,!is.nan(all_merged$Fst_Meadow_Lake))

#MAH Meadow
#MER Lake
#PVM Meadow
#SUM Meadow
#STO Lake
#CHR Lake
#RON Meadow
#ROC Lake
#ELF Lake
#NAM Meadow
#MAR Lake
#PAP Meadow

within=vector(mode="numeric",length=length(all_merged$start))
between=vector(mode="numeric",length=length(all_merged$start))
mean_l_v_l=vector(mode="numeric",length=length(all_merged$start))
mean_l_v_m=vector(mode="numeric",length=length(all_merged$start))
mean_m_v_m=vector(mode="numeric",length=length(all_merged$start))
mean_diff=vector(mode="numeric",length=length(all_merged$start))

for (g in 1:length(all_merged$start)) {
  l_v_l=c(all_merged$Fst_MER_STO[g],all_merged$Fst_MER_CHR[g],all_merged$Fst_MER_ROC[g],
          all_merged$Fst_MER_ELF[g],all_merged$Fst_MER_MAR[g],all_merged$Fst_STO_CHR[g],
          all_merged$Fst_STO_ROC[g],all_merged$Fst_STO_ELF[g],all_merged$Fst_STO_MAR[g],
          all_merged$Fst_CHR_ROC[g],all_merged$Fst_CHR_ELF[g],all_merged$Fst_CHR_MAR[g],
          all_merged$Fst_ROC_ELF[g],all_merged$Fst_ROC_MAR[g],all_merged$Fst_ELF_MAR[g])
  m_v_m=c(all_merged$Fst_MAH_PVM[g],all_merged$Fst_MAH_SUM[g],all_merged$Fst_MAH_RON[g],
          all_merged$Fst_MAH_NAM[g],all_merged$Fst_MAH_PAP[g],all_merged$Fst_PVM_SUM[g],
          all_merged$Fst_PVM_RON[g],all_merged$Fst_PVM_NAM[g],all_merged$Fst_PVM_PAP[g],
          all_merged$Fst_SUM_RON[g],all_merged$Fst_SUM_NAM[g],all_merged$Fst_SUM_PAP[g],
          all_merged$Fst_RON_NAM[g],all_merged$Fst_RON_PAP[g],all_merged$Fst_NAM_PAP[g])
  l_v_m=c(all_merged$Fst_MAH_MER[g],all_merged$Fst_MAH_STO[g],all_merged$Fst_MAH_CHR[g],
          all_merged$Fst_MAH_ROC[g],all_merged$Fst_MAH_ELF[g],all_merged$Fst_MAH_MAR[g],
          all_merged$Fst_MER_PVM[g],all_merged$Fst_MER_SUM[g],all_merged$Fst_MER_RON[g],
          all_merged$Fst_MER_NAM[g],all_merged$Fst_MER_PAP[g],all_merged$Fst_PVM_STO[g],
          all_merged$Fst_PVM_CHR[g],all_merged$Fst_PVM_ROC[g],all_merged$Fst_PVM_ELF[g],
          all_merged$Fst_PVM_MAR[g],all_merged$Fst_SUM_STO[g],all_merged$Fst_SUM_CHR[g],
          all_merged$Fst_SUM_ROC[g],all_merged$Fst_SUM_ELF[g],all_merged$Fst_SUM_MAR[g],
          all_merged$Fst_STO_RON[g],all_merged$Fst_STO_NAM[g],all_merged$Fst_STO_PAP[g],
          all_merged$Fst_CHR_RON[g],all_merged$Fst_CHR_NAM[g],all_merged$Fst_CHR_PAP[g],
          all_merged$Fst_RON_ROC[g],all_merged$Fst_RON_ELF[g],all_merged$Fst_RON_MAR[g],
          all_merged$Fst_ROC_NAM[g],all_merged$Fst_ROC_PAP[g],all_merged$Fst_ELF_NAM[g],
          all_merged$Fst_ELF_PAP[g],all_merged$Fst_NAM_MAR[g],all_merged$Fst_MAR_PAP[g])
  
  mean_l_v_l[g]=mean(l_v_l,na.rm=T)
  mean_m_v_m[g]=mean(m_v_m,na.rm=T)
  mean_l_v_m[g]=mean(l_v_m,na.rm=T)
  
  #combine l_v_l and m_v_m to do ttest below
  combined=c(l_v_l,m_v_m)
  #skip 
  if (all(is.na(l_v_l)) || all(is.na(m_v_m)) || all(is.na(l_v_m)) || sum(!is.na(l_v_l))==1 || sum(!is.na(l_v_m))==1 || sum(!is.na(m_v_m))==1) { 
    skip=all_merged$gene[g]
    print(skip)
    between[g]=NaN
    within[g]=NaN
    mean_diff[g]=NaN
  } else {
    x=t.test(l_v_l,m_v_m)
    within[g]=x$p.value
    y=t.test(combined,l_v_m)
    between[g]=y$p.value
    mean_diff[g]=mean(l_v_m,na.rm=T)-mean(combined,na.rm=T)    
    #fit <- lm(MisSynProp~EcotypeComp, data = dat)
    #slope <- coef(fit)[2]
    #pVal <- anova(fit)$'Pr(>F)'[1]
    
  }
  
}

all_merged$within=within
all_merged$between=between
all_merged$l_v_l_mean=mean_l_v_l
all_merged$l_v_m_mean=mean_l_v_m
all_merged$m_v_m_mean=mean_m_v_m
all_merged$mean_diff=mean_diff


fst_gene_list2=all_merged$INFO.x[all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean & all_merged$mean_diff>0]
unique(fst_gene_list2[!is.na(fst_gene_list2)])

#top_fst=subset(all_merged,all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean)
top_fst2=all_merged[,c(1:4,158:165)]
write.csv(top_fst2,file="per_site_R_outputs/fst_within_v_between_ALL_POPS.csv",quote=F,row.names = F)



#repeat, but with DXY
all_merged=merge(all_fst,fst,by="start",all=T)

#all_merged=subset(all_merged,!is.nan(all_merged$dxy_Meadow_Lake))

within=vector(mode="numeric",length=length(all_merged$start))
between=vector(mode="numeric",length=length(all_merged$start))
mean_l_v_l=vector(mode="numeric",length=length(all_merged$start))
mean_l_v_m=vector(mode="numeric",length=length(all_merged$start))
mean_m_v_m=vector(mode="numeric",length=length(all_merged$start))
mean_diff=vector(mode="numeric",length=length(all_merged$start))

for (g in 1:length(all_merged$start)) {
  l_v_l=c(all_merged$dxy_MER_STO[g],all_merged$dxy_MER_CHR[g],all_merged$dxy_MER_ROC[g],
          all_merged$dxy_MER_ELF[g],all_merged$dxy_MER_MAR[g],all_merged$dxy_STO_CHR[g],
          all_merged$dxy_STO_ROC[g],all_merged$dxy_STO_ELF[g],all_merged$dxy_STO_MAR[g],
          all_merged$dxy_CHR_ROC[g],all_merged$dxy_CHR_ELF[g],all_merged$dxy_CHR_MAR[g],
          all_merged$dxy_ROC_ELF[g],all_merged$dxy_ROC_MAR[g],all_merged$dxy_ELF_MAR[g])
  m_v_m=c(all_merged$dxy_MAH_PVM[g],all_merged$dxy_MAH_SUM[g],all_merged$dxy_MAH_RON[g],
          all_merged$dxy_MAH_NAM[g],all_merged$dxy_MAH_PAP[g],all_merged$dxy_PVM_SUM[g],
          all_merged$dxy_PVM_RON[g],all_merged$dxy_PVM_NAM[g],all_merged$dxy_PVM_PAP[g],
          all_merged$dxy_SUM_RON[g],all_merged$dxy_SUM_NAM[g],all_merged$dxy_SUM_PAP[g],
          all_merged$dxy_RON_NAM[g],all_merged$dxy_RON_PAP[g],all_merged$dxy_NAM_PAP[g])
  l_v_m=c(all_merged$dxy_MAH_MER[g],all_merged$dxy_MAH_STO[g],all_merged$dxy_MAH_CHR[g],
          all_merged$dxy_MAH_ROC[g],all_merged$dxy_MAH_ELF[g],all_merged$dxy_MAH_MAR[g],
          all_merged$dxy_MER_PVM[g],all_merged$dxy_MER_SUM[g],all_merged$dxy_MER_RON[g],
          all_merged$dxy_MER_NAM[g],all_merged$dxy_MER_PAP[g],all_merged$dxy_PVM_STO[g],
          all_merged$dxy_PVM_CHR[g],all_merged$dxy_PVM_ROC[g],all_merged$dxy_PVM_ELF[g],
          all_merged$dxy_PVM_MAR[g],all_merged$dxy_SUM_STO[g],all_merged$dxy_SUM_CHR[g],
          all_merged$dxy_SUM_ROC[g],all_merged$dxy_SUM_ELF[g],all_merged$dxy_SUM_MAR[g],
          all_merged$dxy_STO_RON[g],all_merged$dxy_STO_NAM[g],all_merged$dxy_STO_PAP[g],
          all_merged$dxy_CHR_RON[g],all_merged$dxy_CHR_NAM[g],all_merged$dxy_CHR_PAP[g],
          all_merged$dxy_RON_ROC[g],all_merged$dxy_RON_ELF[g],all_merged$dxy_RON_MAR[g],
          all_merged$dxy_ROC_NAM[g],all_merged$dxy_ROC_PAP[g],all_merged$dxy_ELF_NAM[g],
          all_merged$dxy_ELF_PAP[g],all_merged$dxy_NAM_MAR[g],all_merged$dxy_MAR_PAP[g])
  
  mean_l_v_l[g]=mean(l_v_l,na.rm=T)
  mean_m_v_m[g]=mean(m_v_m,na.rm=T)
  mean_l_v_m[g]=mean(l_v_m,na.rm=T)
  
  #combine l_v_l and m_v_m to do ttest below
  combined=c(l_v_l,m_v_m)
  
  #skip 
  if (all(is.na(l_v_l)) || all(is.na(m_v_m)) || all(is.na(l_v_m)) || sum(!is.na(l_v_l))==1 || sum(!is.na(l_v_m))==1 || sum(!is.na(m_v_m))==1 || sum(!is.na(unique(l_v_l)))==1 || sum(!is.na(unique(l_v_m)))==1 || sum(!is.na(unique(m_v_m)))==1) { 
    skip=all_merged$gene[g]
    print(skip)
    between[g]=NaN
    within[g]=NaN
    mean_diff[g]=NaN
  } else {
    x=t.test(l_v_l,m_v_m)
    within[g]=x$p.value
    y=t.test(combined,l_v_m)
    between[g]=y$p.value
    mean_diff[g]=mean(l_v_m,na.rm=T)-mean(combined,na.rm=T)      
    #fit <- lm(MisSynProp~EcotypeComp, data = dat)
    #slope <- coef(fit)[2]
    #pVal <- anova(fit)$'Pr(>F)'[1]
    
  }
  
}

all_merged$within=within
all_merged$between=between
all_merged$l_v_l_mean=mean_l_v_l
all_merged$l_v_m_mean=mean_l_v_m
all_merged$m_v_m_mean=mean_m_v_m
all_merged$mean_diff=mean_diff

dxy_gene_list2=all_merged$INFO.x[all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean & all_merged$mean_diff>0]
unique(dxy_gene_list2[!is.na(dxy_gene_list2)])

#top_dxy=subset(all_merged,all_merged$between<0.05 & all_merged$l_v_m_mean>all_merged$l_v_l_mean & all_merged$l_v_m_mean>all_merged$m_v_m_mean)
top_dxy2=all_merged[,c(1:4,157,160:165)]
write.csv(top_dxy2,file="per_site_R_outputs/dxy_within_v_between_ALL_POPS.csv",quote=F,row.names = F)

#intersect fst/dxy lists
intersect(fst_gene_list2[!is.na(fst_gene_list2)],dxy_gene_list2[!is.na(dxy_gene_list2)])


#overall intersection
#fst
f=intersect(fst_gene_list2[!is.na(fst_gene_list2)],fst_gene_list[!is.na(fst_gene_list)])
f

#dxy
d=intersect(dxy_gene_list[!is.na(dxy_gene_list)],dxy_gene_list2[!is.na(dxy_gene_list2)])
d

#all sets
intersect(f,d)
#intersect(fst_gene_list2[!is.na(fst_gene_list2)],fst_gene_list[!is.na(fst_gene_list)],dxy_gene_list[!is.na(dxy_gene_list)],dxy_gene_list2[!is.na(dxy_gene_list2)])


# per site FST allele frequencies ------------------------------------------------------------

#setwd("~/Desktop/garter_snakes/missense_results/new/missense/")
#fst=read.csv(file="per_site_shell_outputs/missense-sites.meadow-v-lake.csv",header=T)
#fst=read.table(file="missense-sites.meadow-v-lake.csv",sep="\t",header=T)
#genes=read.csv("pop_gen_summary_lake_v_meadow.csv")
#meadow=read.table(file="per_site_shell_outputs/meadow.missense.frq.count",header=T,row.names=NULL)
#lake=read.table(file="per_site_shell_outputs/lake.missense.frq.count",header=T,row.names=NULL)

#frq.count=merge(lake,meadow,by="CHROM",all=T)

#final=merge(fst,frq.count,by.x="start",by.y="CHROM", all=T)
#final$sample_size=final$N_ALLELES.x+final$N_ALLELES.y
#hist(final$sample_size[final$sample_size<=50])

#colnames(final)=c("start","CHROM","end","mid","sites","pi_Meadow","pi_Lake","dxy_Meadow_Lake","Fst_Meadow_Lake","Gene",
#                  "CHROM.lake","N_ALLELES.lake","N_CHR.lake",
#                  "ALLELE1.COUNT.lake","ALLELE2.COUNT.lake","CHROM.meadow",
#                  "N_ALLELES.meadow","N_CHR.meadow","ALLELE1.COUNT.meadow","ALLELE2.COUNT.meadow","Sample Size")


#write.csv(final[,c(2,1,3,6:10,12:15,17:21)], file="per_site_R_outputs/per_site_fst_and_allele_counts.csv",row.names=F)


#add in allele freqs for ALL pops
MER=read.table(file="per_site_shell_outputs/MER.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
MER=MER[-1,]

PAP=read.table(file="per_site_shell_outputs/PAP.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
PAP=PAP[-1,]

ELF=read.table(file="per_site_shell_outputs/ELF.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
ELF=ELF[-1,]

MAR=read.table(file="per_site_shell_outputs/MAR.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
MAR=MAR[-1,]

NAM=read.table(file="per_site_shell_outputs/NAM.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
NAM=NAM[-1,]

ROC=read.table(file="per_site_shell_outputs/ROC.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
ROC=ROC[-1,]

RON=read.table(file="per_site_shell_outputs/RON.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
RON=RON[-1,]

CHR=read.table(file="per_site_shell_outputs/CHR.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
CHR=CHR[-1,]

PVM=read.table(file="per_site_shell_outputs/PVM.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
PVM=PVM[-1,]

STO=read.table(file="per_site_shell_outputs/STO.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
STO=STO[-1,]

SUM=read.table(file="per_site_shell_outputs/SUM.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
SUM=SUM[-1,]

MAH=read.table(file="per_site_shell_outputs/MAH.frq.count",header=F,row.names=NULL,col.names=c("CHROM","POS","N_ALLELES","N_CHR","ALLELE1.COUNT","ALLELE2.COUNT","ALLELE3.COUNT"),fill=T)
MAH=MAH[-1,]

step1=merge(fst,MER,by.x="start",by.y="POS", all=T)
step2=merge(step1,PAP,by.x="start",by.y="POS", all=T)
step3=merge(step2,ELF,by.x="start",by.y="POS", all=T)
step4=merge(step3,MAR,by.x="start",by.y="POS", all=T)
step5=merge(step4,NAM,by.x="start",by.y="POS", all=T)
step6=merge(step5,ROC,by.x="start",by.y="POS", all=T)
step7=merge(step6,RON,by.x="start",by.y="POS", all=T)
step8=merge(step7,CHR,by.x="start",by.y="POS", all=T)
step9=merge(step8,PVM,by.x="start",by.y="POS", all=T)
step10=merge(step9,STO,by.x="start",by.y="POS", all=T)
step11=merge(step10,SUM,by.x="start",by.y="POS", all=T)
step12=merge(step11,MAH,by.x="start",by.y="POS", all=T)

colnames(step12)=c("start","CHROM","end","mid","sites","pi_Meadow","pi_Lake","dxy_Meadow_Lake","Fst_Meadow_Lake","Gene",
                   "CHROM.MER","N_ALLELES.MER","N_CHR.MER","ALLELE1.COUNT.MER","ALLELE2.COUNT.MER","ALLELE3.COUNT.MER",
                   "CHROM.PAP","N_ALLELES.PAP","N_CHR.PAP","ALLELE1.COUNT.PAP","ALLELE2.COUNT.PAP","ALLELE3.COUNT.PAP",
                   "CHROM.ELF","N_ALLELES.ELF","N_CHR.ELF","ALLELE1.COUNT.ELF","ALLELE2.COUNT.ELF","ALLELE3.COUNT.ELF",
                   "CHROM.MAR","N_ALLELES.MAR","N_CHR.MAR","ALLELE1.COUNT.MAR","ALLELE2.COUNT.MAR","ALLELE3.COUNT.MAR",
                   "CHROM.NAM","N_ALLELES.NAM","N_CHR.NAM","ALLELE1.COUNT.NAM","ALLELE2.COUNT.NAM","ALLELE3.COUNT.NAM",
                   "CHROM.ROC","N_ALLELES.ROC","N_CHR.ROC","ALLELE1.COUNT.ROC","ALLELE2.COUNT.ROC","ALLELE3.COUNT.ROC",
                   "CHROM.RON","N_ALLELES.RON","N_CHR.RON","ALLELE1.COUNT.RON","ALLELE2.COUNT.RON","ALLELE3.COUNT.RON",
                   "CHROM.CHR","N_ALLELES.CHR","N_CHR.CHR","ALLELE1.COUNT.CHR","ALLELE2.COUNT.CHR","ALLELE3.COUNT.CHR",
                   "CHROM.PVM","N_ALLELES.PVM","N_CHR.PVM","ALLELE1.COUNT.PVM","ALLELE2.COUNT.PVM","ALLELE3.COUNT.PVM",
                   "CHROM.STO","N_ALLELES.STO","N_CHR.STO","ALLELE1.COUNT.STO","ALLELE2.COUNT.STO","ALLELE3.COUNT.STO",
                   "CHROM.SUM","N_ALLELES.SUM","N_CHR.SUM","ALLELE1.COUNT.SUM","ALLELE2.COUNT.SUM","ALLELE3.COUNT.SUM",
                   "CHROM.MAH","N_ALLELES.MAH","N_CHR.MAH","ALLELE1.COUNT.MAH","ALLELE2.COUNT.MAH","ALLELE3.COUNT.MAH")

step12=step12[-1,]

write.csv(step12[,c(2,1,3,6:10,12:16,18:22,24:28,30:34,36:40,42:46,48:52,54:58,60:64,66:70,72:76,78:82)], file="per_site_R_outputs/per_site_fst_and_allele_counts_ALL_POPS.csv",row.names=F)


# upset plots -------------------------------------------------------------

library(UpSetR)

AllsitesExpression <- c(lake = 276, meadow = 362, `lake&meadow` = 7118)
missenseExpression <- c(lake = 67, meadow = 87, `lake&meadow` = 1321)

png(file="per_gene_R_outputs/allsites.upset.png",res=150,height=4,width=5,units="in",pointsize=12)
upset(fromExpression(AllsitesExpression), order.by = "freq",mainbar.y.label = "SNP Intersections All Sites",sets.x.label = "SNPs per population")
dev.off()


png(file="per_site_R_outputs/missense.upset.png",res=150,height=4,width=5,units="in",pointsize=12)
upset(fromExpression(missenseExpression), order.by = "freq",mainbar.y.label = "SNP Intersections Missense variants ONLY",sets.x.label = "SNPs per population")
dev.off()

