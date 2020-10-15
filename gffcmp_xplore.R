
library(ggplot2)
library(plotly)
# load a tmap file
tmap<-read.table('data/polyA/cmpR3646.R3646_C28MYACXX_C3VHJACXX.accepted_hits.bam.gtf.tmap', 
               header=T, sep='\t')
#inspect the FPKM and TPM distributions using histograms
ggplotly(ggplot(tmap, aes(x=log2(FPKM)))+geom_histogram(binwidth=0.08))
#to see a nice density line (geom_density() smoothes a histogram):
ggplotly(ggplot(tmap, aes(x=log2(TPM)))+
      geom_histogram(binwidth=0.08, mapping=aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666"))
# note that the geom_density() does not show up unless we use the 
# mapping=aes(y=..density..) option (!?)
# equivalently, we force the density mapping:
ggplotly(ggplot(tmap, aes(x=log2(TPM)))+
           geom_histogram(binwidth=0.08)+geom_density(aes(y=0.08*..count..)))
# filter the data to only multi-exon transcripts and only a few columns of interest
tmmap=fltmap<-subset(tmap, FPKM>1.0 & num_exons>1, 
                    select=c(ref_id,class_code,qry_id, num_exons, FPKM, TPM, len))
#draw a bar chart:
ggplot(tmmap)+geom_bar(mapping=aes(class_code))
#which is the same with
ggplot(tmmap)+geom_bar(aes(x=class_code))
#this is the same with:
ggplot(tmmap, aes(class_code))+geom_bar()
### 1) one way to add labels with actual counts above each bar:
# is to use the data from the bar plot already
g <- ggplot(tmmap)+geom_bar(aes(x=class_code))
g_dt <- layer_data(g) #get layer data from the graph we just built
#then use annotate():
g+annotate(geom="text",label=g_dt$count, x=g_dt$x, y=g_dt$ymax+340)
### 2) another way to add the count labels is to call stat_count() 
### and use the ..count.. internals:
ggplot(tmmap)+geom_bar()+aes(class_code)+
  stat_count(aes(label=..count..), vjust=-0.5, geom="text", position="identity")
# the above is the same with:
ggplot(tmmap, aes(class_code))+geom_bar()+
  stat_count(aes(label=..count..), vjust=-0.5, geom="text", position="identity")
# which is also the same with:
ggplot(tmmap)+aes(class_code) + geom_bar() + 
  geom_text(aes(label=..count..), vjust=-0.5, stat="count", position="identity")
