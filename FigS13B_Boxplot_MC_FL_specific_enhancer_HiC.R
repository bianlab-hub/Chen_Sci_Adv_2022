library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#######MC-specific enhancers Hi-C interaction in MC and FL####################
FL<-read.csv(file = 'MC_specific_enhancers_in_FL_HiC_filtered.bed',sep = '\t',header = F)
FL$Group<-c('FL')
colnames(FL)<-c('Enh_chr','Enh_start','Enh_end','Target_Gene','HiC_interaction','Group')
FL <- FL[!duplicated(FL$Enh_start),]
MC<-read.csv(file = 'MC_specific_enhancers_in_MC_HiC_filtered.bed',sep = '\t',header = T)
MC$Group<-c('MC')
colnames(MC)<-c('Enh_chr','Enh_start','Enh_end','Target_Gene','HiC_interaction','Group')
MC <- MC[!duplicated(MC$Enh_start),]
Boxplot<-rbind(MC,FL)
cols <- c("Blue", "Red")
pdf('HiC_MC-specific_enhancers_in_MC_FL.pdf', width=3, height=4)
ggplot(Boxplot,aes(x=factor(Group,levels = c('MC','FL')),y=HiC_interaction,color=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(0,0.02))+labs(y = c('HiC_interaction'))+
  stat_compare_means(comparisons = list(c('MC','FL')),
                     method="t.test")
dev.off()
ttest<-t.test(Boxplot$HiC_interaction~Boxplot$Group)
ttest$p.value
write.table(ttest$p.value,file = 'HiC_MC-specific_enhancers_in_MC_FL_ttest.txt',col.names = F,row.names = F,quote = F)
#######FL-specific enhancers Hi-C interaction in MC and FL####################
FL<-read.csv(file = 'FL_specific_enhancers_in_FL_HiC_filtered.bed',sep = '\t',header = F)
FL$Group<-c('FL')
colnames(FL)<-c('Enh_chr','Enh_start','Enh_end','Target_Gene','HiC_interaction','Group')
FL <- FL[!duplicated(FL$Enh_start),]
MC<-read.csv(file = 'FL_specific_enhancers_in_MC_HiC_filtered.bed',sep = '\t',header = F)
MC$Group<-c('MC')
colnames(MC)<-c('Enh_chr','Enh_start','Enh_end','Target_Gene','HiC_interaction','Group')
MC <- MC[!duplicated(MC$Enh_start),]
Boxplot<-rbind(MC,FL)
cols <- c("Blue", "Red")
pdf('HiC_FL-specific_enhancers_in_MC_FL.pdf', width=3, height=4)
ggplot(Boxplot,aes(x=factor(Group,levels = c('MC','FL')),y=HiC_interaction,color=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(0,0.02))+labs(y = c('HiC_interaction'))+
  stat_compare_means(comparisons = list(c('MC','FL')),
                     method="t.test")
dev.off()
ttest<-t.test(Boxplot$HiC_interaction~Boxplot$Group)
ttest$p.value
write.table(ttest$p.value,file = 'HiC_FL-specific_enhancers_in_MC_FL_ttest.txt',col.names = F,row.names = F,quote = F)