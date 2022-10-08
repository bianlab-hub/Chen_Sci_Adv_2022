fl <- read.csv('forelimb_chr11_112.75_112.8_chr11_111_114_balanced.tsv', sep='\t', header=TRUE)

mc <- read.csv('mc_chr11_112.75_112.8_chr11_111_114_balanced.tsv', sep='\t', header=TRUE)

x_coordinate <- seq(from=112755000, to=112795000, by =10000)
y_coordinate <- seq(from=111005000, to=113995000, by =10000)

colnames(fl)<- y_coordinate
colnames(fl)
rownames(fl) <- x_coordinate

colnames(mc)<- y_coordinate
rownames(mc) <- x_coordinate


# pseudo 4C using 10kb window of 112.77-78 as viewpoint
fl_sox9 <- fl[3,]/ sum(fl[3,], na.rm = TRUE)
mc_sox9 <- mc[3,]/ sum(mc[3,], na.rm = TRUE)
mc_diff_fl <- mc_sox9 - fl_sox9



pdf('sox9_pseudo4c_FL.pdf', width=10, height=5)
plot(seq(from=111005000, to=113995000, by =10000), fl_sox9, col='white', ylab='Interaction Frequency',xlab='Chr11 Coordinate',
     main='Sox9 Pseudo 4C', type='l',cex.lab=0.8,cex.main=0.8,cex.axis=0.8,xlim=c(111500000, 113500000), ylim=c(0, 0.03))
lines(seq(from=111005000, to=113995000, by =10000), fl_sox9 , lwd=2, col='blue')
legend('topright', legend=c('FL') ,col=c('blue') , lwd=2 , bty='n', cex=2)
## PRS region
abline(v=c(111555526,111833044),col='black',lty=2)
## Proximal TBC
abline(v=c(112437850,112722603),col='grey',lty=2)
dev.off()

pdf('sox9_pseudo4c_MC.pdf', width=10, height=5)
plot(seq(from=111005000, to=113995000, by =10000), mc_sox9, col='white', ylab='Interaction Frequency',xlab='Chr11 Coordinate',
     main='Sox9 Pseudo 4C', type='l',cex.lab=0.8,cex.main=0.8,cex.axis=0.8,xlim=c(111500000, 113500000), ylim=c(0, 0.03))
lines(seq(from=111005000, to=113995000, by =10000), fl_sox9 , lwd=2, col='blue')
legend('topright', legend=c('MC') ,col=c('red') , lwd=2 , bty='n', cex=2)
## PRS region
abline(v=c(111555526,111833044),col='black',lty=2)
## Proximal TBC
abline(v=c(112437850,112722603),col='grey',lty=2)
dev.off()
