library(ggplot2)
library(plotly)
library(dplyr)
library(cowplot)

# load the aligner data results
t <- '2'
adf<-read.table('txnoise/results.head.tab', 
                  header=T, sep='\t')
for (t in 0:2) {
 for (s in 0:9) {
   bdf<-read.table(paste0('txnoise/t',toString(t),'_s',toString(s),'/results.tab'), 
                header=F, sep='\t')
   colnames(bdf) <- colnames(adf)
   adf <- rbind(adf, bdf)
  }
}

at <- subset(adf, grepl('_s\\d$', adf$sample))
ats <- subset(adf, grepl('\\.s$', adf$sample))
#% percent of correct spliced alignments 
ggplot(data=subset(at, aligner %in% c('hisat2','star','hisat2nodta','hisat21')), aes(x=as.numeric(factor(sample)), y=num_match_spl/num_sim_spl, group=aligner, color=aligner)) +
  geom_line(size=1) + geom_point(size=1.5) + xlab("samples") + ylab("% correct spliced alignments")+ggtitle("Percentage of correct spliced alignments")+theme(plot.title = element_text(hjust = 0.5))
# ggplot(data=subset(at, aligner %in% c('hisat2','star','hisat2nodta','hisat21')),
# aes(x=as.numeric(factor(sample)), y=num_match_spl/num_sim_spl, 
# group=aligner, color=aligner)) + geom_line(size=1.2) + geom_point(size=2) + xlab("samples") + ylab("% correct spliced alignments")
# % correct alignments overall (spliced and unspliced)
ggplot(data=at, aes(x=as.numeric(factor(sample)), y=num_match/num_sim, group=aligner, color=aligner)) +
  geom_line(size=1) + geom_point(size=1.5) + xlab("samples") + ylab("% correct alignments")+ggtitle("Percentage of correct alignments")+theme(plot.title = element_text(hjust = 0.5))

cmpdf <- read.table('txnoise/t0_s0/gffcmp_summary.tab', 
                    header=T, sep='\t', nrows=1)
cmpdf <- cmpdf[0,]
colnames(cmpdf) <- c('sample', colnames(cmpdf))
for (t in 0:2) {
  for (s in 0:9) {
    smp <- paste0('t',toString(t),'_s',toString(s))
    bdf<-read.table(paste0('txnoise/',smp,'/gffcmp_summary.tab'), 
                    header=T, sep='\t')
    #colnames(bdf) <- colnames(adf)
    cmpdf <- rbind(cmpdf, bdf)
  }
}

ggplot(data=subset(cmpdf, queryFileName=='test.strg'), aes(x=as.numeric(factor(sample)), y=num_match/num_sim, group=aligner, color=aligner)) +
  geom_line(size=1) + geom_point(size=1.5) + xlab("samples") + ylab("% correct alignments")+ggtitle("Percentage of correct alignments")+theme(plot.title = element_text(hjust = 0.5))

  
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



