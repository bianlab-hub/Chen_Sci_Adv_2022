library(dplyr)
library(ggplot2)
MC_insulation <- read.csv('MC_HiC_combined_50kb_insulation', sep='\t', header=T)
MC_insulation_chr11<-dplyr::filter(MC_insulation,chrom == 'chr11')
FL_insulation <- read.csv('FL_HiC_combined_50kb_insulation', sep='\t', header=T)
FL_insulation_chr11<-dplyr::filter(FL_insulation,chrom == 'chr11')

pdf(file = 'Sox9_Kcnj2_TAD_insulation_score_MC_vs_FL.pdf',width = 12,height = 6)
ggplot()+geom_line(data = MC_insulation_chr11,aes(x = start, y = log2_insulation_score_500000,colour = "MC_insulation"),size=0.4)+
  geom_line(data = FL_insulation_chr11,aes(x = start, y = log2_insulation_score_500000,colour = "FL_insulation"),size=0.4)+
  xlab("chr11 Coordinates")+ylab("Insulation Score")+ylim(-1.5,1.5)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  xlim(c(110000000,114000000))+
  scale_colour_manual("",values = c("MC_insulation" = "red","FL_insulation"="blue"))+
  ggtitle("Sox9-Kcnj2 locus")+geom_vline(xintercept = c(111555526,111833044),linetype = 'dashed')
dev.off()
