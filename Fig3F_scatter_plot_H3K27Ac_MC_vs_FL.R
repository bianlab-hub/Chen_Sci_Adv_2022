library(ggplot2)
scatter<- read.csv('H3K27Ac_input_for_scatter_plot.csv',header = T)
cols <- c("Grey", "Blue", "Red")
log10MC<-log10(scatter$Normalized_H3K27Ac_MC)
log10FL<-log10(scatter$Normalized_H3K27Ac_FL)
## Scatter Plot based on Normalized H3K27Ac with smoothing and regression
pdf(file = "H3K27Ac_MC_vs_FL_smoothed_log10_regression_all.pdf",height = 5,width = 6)
ggplot(data = scatter)+
  geom_smooth(aes(x=log10MC,y=log10FL),method = "lm",se = F, color = 'black')+
  scale_color_manual(values = cols)+
  geom_point(aes(x=log10MC,y=log10FL,colour=Region))+
  labs(y = "log10(Normalized H3K27Ac FL)",
       x="log10(Normalized H3K27Ac MC)")
dev.off()

