library(ggplot2)
#library(plotly)
library(dplyr)
library(cowplot)

# load a tmap file
cmpdf<-read.table('new_gffcmp_stats.tab', 
                  header=T, sep='\t')
cmpA<-subset(cmpdf, protocol=="polyA")
cmpZ<-subset(cmpdf, protocol=="riboZ")
#transcript sensitivity 
tSnA<-ggplot(cmpA, aes(fill=aligner, y=tSn, x=qryfile)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=tSn),position=position_dodge(width=0.9),vjust=-0.25)+
  coord_cartesian(ylim = c(min(cmpA$tSn),max(cmpA$tSn)+0.3))+
  labs(x="polyA samples", y="Sensitivity")

tSnZ<-ggplot(cmpZ, aes(fill=aligner, y=tSn, x=qryfile)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=tSn),position=position_dodge(width=0.9),vjust=-0.25)+
  coord_cartesian(ylim = c(min(cmpZ$tSn),max(cmpZ$tSn)+0.3))+
  labs(x="riboZero samples", y="Sensitivity")
plot_grid(tSnA, tSnZ, nrow=2, ncol=1)

#transcript precision
tPrA<-ggplot(cmpA, aes(fill=aligner, y=tPr, x=qryfile)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=tPr),position=position_dodge(width=0.9),vjust=-0.25)+
  coord_cartesian(ylim = c(min(cmpA$tPr),max(cmpA$tPr)+0.3))+
  labs(x="polyA samples", y="Precision")

tPrZ<-ggplot(cmpZ, aes(fill=aligner, y=tPr, x=qryfile)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label=tPr),position=position_dodge(width=0.9),vjust=-0.25)+
  coord_cartesian(ylim = c(min(cmpZ$tPr),max(cmpZ$tPr)+0.3))+
  labs(x="riboZero samples", y="Precision")
plot_grid(tPrA, tPrZ, nrow=2, ncol=1)



