library(ggplot2)
library(ggpubr)
library(RColorBrewer)

Boxplot<-read.csv('MC_vs_FL_ATAC_Boxplot_input_PRS_TBC.csv', sep=',', header=TRUE)
cols <- c("Blue", "Red")
pdf('ATAC_PRS_vs_TBC_log2.pdf', width=3, height=6)
ggplot(Boxplot,aes(x=Region,y=ATAC_log2,color=Region))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(-1,1))+
  geom_jitter(aes(x=Region,y=ATAC_log2,color=Region))+
  stat_compare_means(comparisons = list(c('PRS_region','Proximal_TBC')                                        ),
                     method="t.test",margin_top = 0.05)
dev.off()

Boxplot<-read.csv('MC_vs_FL_H3K27Ac_Boxplot_input_PRS_TBC.csv', sep=',', header=TRUE)
cols <- c("Blue", "Red")
pdf('H3K27Ac_PRS_vs_TBC_log2.pdf', width=3, height=6)
ggplot(Boxplot,aes(x=Region,y=H3K27Ac_log2,color=Region))+
  geom_boxplot()+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )+
  scale_color_manual(values = cols)+ylim(c(-1,1))+
  geom_jitter(aes(x=Region,y=H3K27Ac_log2,color=Region))+
  stat_compare_means(comparisons = list(c('PRS_region','Proximal_TBC')                                        ),
                     method="t.test",margin_top = 0.05)
dev.off()

