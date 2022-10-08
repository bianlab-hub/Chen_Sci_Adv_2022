library('ggplot2')
library('dplyr')
ABC_score<-read.csv(file = 'MC_vs_FL_ABC_score_contact_corrected_full.csv')
ABC_score<-ABC_score[,c(1,2,3,4,10)]
ABC_score<-dplyr::filter(ABC_score, MC_end < 112722603) # Remove candidate enhancers near Sox9 promoter
colnames(ABC_score)<-c('chr','start','end','Enh','MC_minus_FL')
PRS_enhs<-dplyr::filter(ABC_score, start > 111555526 & end < 111833044 )
TBC_enhs<-dplyr::filter(ABC_score, start > 112437850 & end < 112722603 )

##########################Set threshold for ABC score################
# MC-specific enhancers###################
##P value
for (i in 1:20){
  MC_specific<-dplyr::filter(ABC_score, MC_minus_FL > (i/2000) )
  MC_specific_in_PRS<-dplyr::filter(ABC_score, MC_minus_FL > (i/2000) & start > 111555526 & end < 111833044)
  Chi_table<-matrix(c(nrow(MC_specific_in_PRS),
                      nrow(PRS_enhs)-nrow(MC_specific_in_PRS),
                      nrow(MC_specific)-nrow(MC_specific_in_PRS),
                      nrow(ABC_score)-nrow(PRS_enhs)-(nrow(MC_specific)-nrow(MC_specific_in_PRS))),
                      nrow = 2,ncol=2)
  test<-mcnemar.test(Chi_table)
  print(i/2000)
  print(Chi_table)
  print(test$p.value)
}
##Predicted MC-specific Enhancers fall within PRS region
for (i in 1:20){
  MC_specific<-dplyr::filter(ABC_score, MC_minus_FL > (i/2000) )
  MC_specific_in_PRS<-dplyr::filter(ABC_score, MC_minus_FL > (i/2000) & start > 111555526 & end < 111833044)
  Chi_table<-matrix(c(nrow(MC_specific_in_PRS),
                      nrow(PRS_enhs)-nrow(MC_specific_in_PRS),
                      nrow(MC_specific)-nrow(MC_specific_in_PRS),
                      nrow(ABC_score)-nrow(PRS_enhs)-(nrow(MC_specific)-nrow(MC_specific_in_PRS))),
                    nrow = 2,ncol=2)
  test<-mcnemar.test(Chi_table)
  PRS_rate<-nrow(MC_specific_in_PRS)/nrow(PRS_enhs)
  print(PRS_rate)
}
##Predicted MC-specific Enhancers fall within region other than PRS region
for (i in 1:20){
  MC_specific<-dplyr::filter(ABC_score, MC_minus_FL > (i/2000) )
  MC_specific_in_PRS<-dplyr::filter(ABC_score, MC_minus_FL > (i/2000) & start > 111555526 & end < 111833044)
  Chi_table<-matrix(c(nrow(MC_specific_in_PRS),
                      nrow(PRS_enhs)-nrow(MC_specific_in_PRS),
                      nrow(MC_specific)-nrow(MC_specific_in_PRS),
                      nrow(ABC_score)-nrow(PRS_enhs)-(nrow(MC_specific)-nrow(MC_specific_in_PRS))),
                    nrow = 2,ncol=2)
  test<-mcnemar.test(Chi_table)
  non_PRS_rate<-(nrow(MC_specific)-nrow(MC_specific_in_PRS))/(nrow(ABC_score)-nrow(PRS_enhs))
  print(non_PRS_rate)
}
ABC_thresholding_MC<-read.csv(file = 'ABC_thresholding_MC.csv',sep = ',')
ggplot(data = ABC_thresholding_MC, 
       mapping = aes(x = Threshold, y = MC_specific_in_PRS)) +geom_line()+geom_point()

pdf(file = 'ABC_thresholding_MC_specific_enhancers.pdf',width = 16,height = 8)
ggplot()+geom_line(data = ABC_thresholding_MC,aes(x = Threshold,y = MC_specific_in_PRS,colour = "PRS_region"),size=0.5)+
  geom_line(data = ABC_thresholding_MC,aes(x = Threshold,y = MC_specific_not_in_PRS,colour = "non-PRS_region"),size=0.5)+
  scale_colour_manual("",values = c("PRS_region"='red',"non-PRS_region"='grey'))+
  geom_point()+
  xlab("ABC threshold")+ylab("Rate of enhancers above ABC threshold")+
  scale_x_continuous(breaks = ABC_thresholding_MC$Threshold)+
  geom_text(data = ABC_thresholding_MC,aes(x = Threshold,
                                           y = MC_specific_in_PRS,
                                           label=P_value,
                                           vjust=-0.5,
                                           hjust=0.5,angle=30),
            size = 2.5)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = c('none')
  )
dev.off()
