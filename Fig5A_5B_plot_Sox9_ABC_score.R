library('ggplot2')
ABC_Barplot <- read.csv('ABC_Score_Barplot.csv', sep=',', header=TRUE)
pdf(file = 'MC_vs_FL_ABC_score_barplot.pdf',height = 5,width = 10)
ggplot(data= ABC_Barplot, aes(x=Enh_Mid, y=ABC_Score, fill=Group)) +
  geom_bar(stat="identity")+ylim(-0.04,0.04)+
  scale_fill_manual(values = c("FL" = "blue","MC"="red"))+
  geom_vline(xintercept = 111555526,linetype = "dashed",size = 0.2)+
  geom_vline(xintercept = 111833044,linetype = "dashed",size = 0.2)+
  geom_vline(xintercept = 112437850,linetype = "dashed",size = 0.2)+
  geom_vline(xintercept = 112722603,linetype = "dashed",size = 0.2)
dev.off()

ABC_Barplot_diff <- read.csv('ABC_Score_diff_Barplot.csv', sep=',', header=TRUE)
pdf(file = 'MC_vs_FL_ABC_score_diff_barplot.pdf',height = 5,width = 10)
ggplot(data= ABC_Barplot_diff, aes(x=Enh_Mid, y=MC_vs_FL)) +
  geom_bar(stat="identity",fill="black")+
  geom_vline(xintercept = 111555526,linetype = "dashed",size = 0.2)+
  geom_vline(xintercept = 111833044,linetype = "dashed",size = 0.2)+
  geom_vline(xintercept = 112437850,linetype = "dashed",size = 0.2)+
  geom_vline(xintercept = 112722603,linetype = "dashed",size = 0.2)+
  geom_hline(yintercept = c(0.0035,-0.0035),linetype='dashed',size=0.2)
dev.off()

pdf(file = 'MC_vs_FL_ABC_score_diff_barplot_Enh38_Enh191.pdf',height = 5,width = 10)
ggplot(data= ABC_Barplot_diff, aes(x=Enh_Mid, y=MC_vs_FL)) +
  geom_bar(stat="identity",fill="black")+
  geom_vline(xintercept = 111567283,linetype = "dashed",size = 0.1)+
  geom_vline(xintercept = 111567783,linetype = "dashed",size = 0.1)+
  geom_vline(xintercept = 112532081,linetype = "dashed",size = 0.1)+
  geom_vline(xintercept = 112532581,linetype = "dashed",size = 0.1)+
  geom_hline(yintercept = c(0.0035,-0.0035),linetype='dashed',size=0.2)
dev.off()
