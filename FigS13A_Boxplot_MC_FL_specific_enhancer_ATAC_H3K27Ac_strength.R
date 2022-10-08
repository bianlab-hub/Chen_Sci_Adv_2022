library(ggplot2)
library(ggpubr)
library(RColorBrewer)
#ATAC
###################FL-specific enhancers in MC and FL##############################
FL<-read.csv(file = 'FL_ATAC_strength_at_FL_specific_enhancers.bed',sep = '\t',header = F)
FL$Group<-c('FL')
colnames(FL)<-c('Enh_chr','Enh_start','Enh_end','ATAC_strength','Group')
FL <- FL[!duplicated(FL$Enh_start),]
MC<-read.csv(file = 'MC_ATAC_strength_at_FL_specific_enhancers.bed',sep = '\t',header = F)
MC$Group<-c('MC')
colnames(MC)<-c('Enh_chr','Enh_start','Enh_end','ATAC_strength','Group')
MC <- MC[!duplicated(MC$Enh_start),]
Boxplot<-rbind(MC,FL)
cols<-c('Blue','Red')
pdf('ATAC_FL-specific_enhancers_in_MC_FL.pdf', width=3, height=4)
ggplot(Boxplot,aes(x=factor(Group,levels = c('MC','FL')),y=ATAC_strength,color=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(0,25))+labs(y = c('ATAC strength'))+
  stat_compare_means(comparisons = list(c('MC','FL')),
                     method="t.test")
dev.off()
ttest<-t.test(Boxplot$ATAC_strength~Boxplot$Group)
ttest$p.value
write.table(ttest$p.value,file = 'FL_ATAC_strength_at_FL_specific_enhancers_ttest.txt',col.names = F,row.names = F,quote = F)
###################MC-specific enhancers in MC and FL##############################
FL<-read.csv(file = 'FL_ATAC_strength_at_MC_specific_enhancers.bed',sep = '\t',header = F)
FL$Group<-c('FL')
colnames(FL)<-c('Enh_chr','Enh_start','Enh_end','ATAC_strength','Group')
FL <- FL[!duplicated(FL$Enh_start),]
MC<-read.csv(file = 'MC_ATAC_strength_at_MC_specific_enhancers.bed',sep = '\t',header = F)
MC$Group<-c('MC')
colnames(MC)<-c('Enh_chr','Enh_start','Enh_end','ATAC_strength','Group')
MC <- MC[!duplicated(MC$Enh_start),]
Boxplot<-rbind(MC,FL)
cols<-c('Blue','Red')
pdf('ATAC_MC-specific_enhancers_in_MC_FL.pdf', width=3, height=4)
ggplot(Boxplot,aes(x=factor(Group,levels = c('MC','FL')),y=ATAC_strength,color=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(0,25))+labs(y = c('ATAC strength'))+
  stat_compare_means(comparisons = list(c('MC','FL')),
                     method="t.test")
dev.off()
ttest<-t.test(Boxplot$ATAC_strength~Boxplot$Group)
ttest$p.value
write.table(ttest$p.value,file = 'MC_ATAC_strength_at_FL_specific_enhancers_ttest.txt',col.names = F,row.names = F,quote = F)



#H3K27Ac
###################FL-specific enhancers in MC and FL##############################
FL<-read.csv(file = 'FL_H3K27Ac_strength_at_FL_specific_enhancers.bed',sep = '\t',header = F)
FL$Group<-c('FL')
colnames(FL)<-c('Enh_chr','Enh_start','Enh_end','H3K27Ac_strength','Group')
FL <- FL[!duplicated(FL$Enh_start),]
MC<-read.csv(file = 'MC_H3K27Ac_strength_at_FL_specific_enhancers.bed',sep = '\t',header = F)
MC$Group<-c('MC')
colnames(MC)<-c('Enh_chr','Enh_start','Enh_end','H3K27Ac_strength','Group')
MC <- MC[!duplicated(MC$Enh_start),]
Boxplot<-rbind(MC,FL)
cols<-c('Blue','Red')
pdf('H3K27Ac_FL-specific_enhancers_in_MC_FL.pdf', width=3, height=4)
ggplot(Boxplot,aes(x=factor(Group,levels = c('MC','FL')),y=H3K27Ac_strength,color=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(0,100))+labs(y = c('H3K27Ac strength'))+
  stat_compare_means(comparisons = list(c('MC','FL')),
                     method="t.test")
dev.off()
ttest<-t.test(Boxplot$H3K27Ac_strength~Boxplot$Group)
ttest$p.value
write.table(ttest$p.value,file = 'FL_H3K27Ac_strength_at_FL_specific_enhancers_ttest.txt',col.names = F,row.names = F,quote = F)
###################MC-specific enhancers in MC and FL##############################
FL<-read.csv(file = 'FL_H3K27Ac_strength_at_MC_specific_enhancers.bed',sep = '\t',header = F)
FL$Group<-c('FL')
colnames(FL)<-c('Enh_chr','Enh_start','Enh_end','H3K27Ac_strength','Group')
FL <- FL[!duplicated(FL$Enh_start),]
MC<-read.csv(file = 'MC_H3K27Ac_strength_at_MC_specific_enhancers.bed',sep = '\t',header = F)
MC$Group<-c('MC')
colnames(MC)<-c('Enh_chr','Enh_start','Enh_end','H3K27Ac_strength','Group')
MC <- MC[!duplicated(MC$Enh_start),]
Boxplot<-rbind(MC,FL)
cols<-c('Blue','Red')
pdf('H3K27Ac_MC-specific_enhancers_in_MC_FL.pdf', width=3, height=4)
ggplot(Boxplot,aes(x=factor(Group,levels = c('MC','FL')),y=H3K27Ac_strength,color=Group))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(0,100))+labs(y = c('H3K27Ac strength'))+
  stat_compare_means(comparisons = list(c('MC','FL')),
                     method="t.test")
dev.off()
ttest<-t.test(Boxplot$H3K27Ac_strength~Boxplot$Group)
ttest$p.value
write.table(ttest$p.value,file = 'MC_H3K27Ac_strength_at_FL_specific_enhancers_ttest.txt',col.names = F,row.names = F,quote = F)


