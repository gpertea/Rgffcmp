
library(ggplot2)
library(plotly)
# load a tmap file
tmap<-read.table('data/riboZ/cmpR3907.R3907_C3RGLACXX.gtf.tmap', 
               header=T, sep='\t')
gtitle <- "R3907 riboZ"
##inspect the FPKM and TPM distributions using histograms
#ggplotly(ggplot(tmap, aes(x=log2(FPKM)))+geom_histogram(binwidth=0.08))
#TPM plot: let's add a nice density line (geom_density() which smoothes the histogram):
ggplotly(ggplot(tmap, aes(x=log2(TPM)))+
           geom_histogram(binwidth=0.1, mapping=aes(y=..density..))
         +stat_density(color="#4060ff", geom="line")+ggtitle(gtitle))
#or, to see density as transparent overlay:
#ggplotly(ggplot(tmap, aes(x=log2(TPM)))+
#      geom_histogram(binwidth=0.08, mapping=aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666"))
# note that the geom_density() does not show up unless we use the 
# mapping=aes(y=..density..) option (!?)
# equivalently, we force the density mapping:
#ggplotly(ggplot(tmap, aes(x=log2(TPM)))+
#           geom_histogram(binwidth=0.08)+geom_density(color='#FF0000', aes(y=0.08*..count..)))+stat_density(geom="line")

##------------------------------


# filter the data to only multi-exon transcripts and only a few columns of interest
tmmap<-subset(tmap, FPKM>1.0 & num_exons>1, 
                    select=c(ref_id,class_code,qry_id, num_exons, FPKM, TPM, len))
##draw a bar chart:
#ggplot(tmmap)+geom_bar(mapping=aes(class_code))
##which is the same with
#ggplot(tmmap)+geom_bar(aes(x=class_code))
##this is the same with:
#ggplot(tmmap, aes(class_code))+geom_bar()
### 1) one way to add labels with actual counts above each bar:
# is to use the data from the bar plot already
g <- ggplot(tmmap)+geom_bar(aes(x=class_code))
g_dt <- layer_data(g) #get layer data from the graph we just built
#then use annotate():
g+annotate(geom="text",label=g_dt$count, x=g_dt$x, y=g_dt$ymax+340)
### 2) another way to add the count labels is to call stat_count() 
### and use the ..count.. internals:
#ggplot(tmmap)+geom_bar()+aes(class_code)+
#  stat_count(aes(label=..count..), vjust=-0.5, geom="text", position="identity")
## the above is the same with:
#ggplot(tmmap, aes(class_code))+geom_bar()+
#  stat_count(aes(label=..count..), vjust=-0.5, geom="text", position="identity")
## which is also the same with:
#ggplot(tmmap)+aes(class_code) + geom_bar() + 
#  geom_text(aes(label=..count..), vjust=-0.5, stat="count", position="identity")
## to hide the legend, colorize the bars and increase font:
#ggplot(tmmap)+aes(class_code) + geom_bar(aes(fill=factor(class_code))) + 
#  geom_text(aes(label=..count..), vjust=-0.5, stat="count", position="identity")+
#  theme(text = element_text(size=16), legend.position = "none")
##
## Want to discard all class codes with less than 100 entries (u,y,o,i etc.)
## while also changing the bar chart so the bars are sorted by frequency
## Note that levels(factor(tmmap$class_code)) has the default order of the bars
## so we have to change that
## -- Create a new column with the desired ordering as given by table(tmap$class_code)
## sort(table(tmmap$class_code), decreasing=TRUE) shows the desired order!
## take names() of that to get just the list of class_codes in the order we want them, 
## so the new factor column should be built like this
## - we could sort by frequency
#> tmmap$class <- factor(tmmap$class_code, 
#                         levels=names(sort(table(tmmap$class_code), decreasing=TRUE)))
### -- but let's use order or interest: 
codePrio <- c('=','j','k','y','i','u','o','x','m','n','s')
tmmap$class <- factor(tmmap$class_code, codePrio)
## now the class column is of type factor, it shows the same class_code 
## but the factor levels internally are different
## the following should now show the class codes (class entries) in the order of 
## decreasing frequency but also discarded if they are lower than 100
#> levels(tmmap$class)[table(tmmap$class)>=100]
## so now we can subset the data frame using this new class column:
# tmsel <- subset(tmmap, class %in% levels(tmmap$class)[table(tmmap$class)>=100])
tmsel <- subset(tmmap, class %in% levels(tmmap$class))
# note the above is the same with this (the class only matters for sorting!):
#tmsel <- subset(tmmap, class_code %in% levels(factor(tmmap$class_code))[table(tmmap$class_code)>=100])
gp1 <- ggplot(tmsel)+aes(class) + geom_bar(aes(fill=class)) + 
  geom_text(aes(label=..count..), vjust=-0.5, stat="count", position="identity")+
  theme(text = element_text(size=16), legend.position = "none") + ggtitle(gtitle)
#ggplotly(gp1)
gp1

