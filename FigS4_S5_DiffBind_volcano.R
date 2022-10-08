library(magrittr)
library(ggpubr)
library(ggthemes)
deg.data<- read.table("Forelimb_vs_MC_edgeR.txt",header = T,sep = "\t")
head(deg.data)
deg.data$logP <- -log10(deg.data$FDR)
deg.data$Group = "not-significant"
deg.data$Group[which( (deg.data$FDR <0.05) & (deg.data$Fold > 1) )]="Higher strength in MC"
deg.data$Group[which( (deg.data$FDR <0.05) & (deg.data$Fold < -1) )]="Higher strength in FL"
table(deg.data$Group)
deg.data$Label = ""
deg.data <- deg.data[order(deg.data$FDR), ]
up.genes <- head(deg.data$Symbol[which(deg.data$Group == "up-regulated")], 10000)
down.genes <- head(deg.data$Symbol[which(deg.data$Group == "down-regulated")], 10000)
deg.top10.genes <- c(as.character(up.genes), as.character(down.genes))
deg.data$Label[match(deg.top10.genes,deg.data$Symbol)] <- deg.top10.genes
deg.volcano <- ggscatter(deg.data,x = "Fold", y="logP", color="Group",
                         palette = c("#0000ff","#ff0000", "#BBBBBB"),
                         size = 1,
                         label = deg.data$Label,
                         font.label = 8,
                         repel = T,
                         xlab = "log2(FoldChange) MC vs FL",
                         ylab = "-log10(FDR)") + theme_base()+
                        geom_vline(xintercept = c(-1,1), linetype = "dashed")+
  geom_hline(yintercept = 1.3, linetype = "dashed")+xlim(-6.5,6.5)
ggsave(deg.volcano,filename = 'deg.volcano.pdf',width = 20, height = 20, units = "cm")

write.csv(deg.data,file = 'H3K27Ac_whole_genome_DiffBind_results.csv',quote = F,sep = ',')
