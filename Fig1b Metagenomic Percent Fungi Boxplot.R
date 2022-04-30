library("ggplot2"); packageVersion("ggplot2")
library("plyr")

AllWGSpercent <- read.table(file="/Users/tommyauchtung/Desktop/Fig1b.txt", header = TRUE, sep = "\t")

############## FIG 1B PORTION THAT HAS TEDDY ###################
justTEDDY <- subset(AllWGSpercent, Project!="3_HMP" & Project!="NonWestern")
ggplot(justTEDDY, aes(x=Project,y=LogPercent)) +
  geom_point(aes(color=Country),position=position_jitterdodge(dodge.width = 0.7)) +
  geom_boxplot(aes(fill=Country),  fatten = 4, alpha=0.7, width=0.6, position = position_dodge(width = 0.7), outlier.shape=NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  axis.line    = element_line(color='black'),
  legend.title = element_blank(),
  legend.position = "none",
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin=margin(l=20,b=15,t=15,r=40,unit = "pt")) +
  scale_fill_manual(values=c("#002f6c","#FFCD00","#FF0000","lawngreen","#002f6c","#FFCD00","#FF0000","lawngreen")) +
  scale_colour_manual(values=c("#000000","#000000","#000000","#000000","#000000","#000000","#000000","#000000")) +
  scale_y_continuous(limits=c(-6,2), breaks=c(-6,-5,-4,-3,-2,-1,0,1,2)) +
  coord_cartesian(ylim = c(-6,2), expand = FALSE)
#Export as 12 x 20 inch pdf

############## FIG 1B PORTION THAT HAS HMP ###################
HMPpercent <- subset(AllWGSpercent, Study!="TEDDY" & Study!="NonWestern")
ggplot(HMPpercent, aes(x=Project,y=LogPercent)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) +
  geom_boxplot(fatten = 4) +
  geom_jitter(width = 0.2) +
  scale_y_continuous(limits=c(-6,2), breaks=c(-6,-5,-4,-3,-2,-1,0,1,2))+
  coord_cartesian(ylim = c(-6,2), expand = FALSE)
#Export as 5 x 20 inch pdf

############## FIG 1B PORTION THAT HAS NON-WESTERN DATA ###################
NWpercent <- subset(AllWGSpercent, Study!="TEDDY" & Study!="HMP")
ggplot(NWpercent, aes(x=Project,y=LogPercent)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) +
  geom_boxplot(fatten = 4) +
  geom_jitter(width = 0.2) +
  scale_y_continuous(limits=c(-6,2), breaks=c(-6,-5,-4,-3,-2,-1,0,1,2))+
  coord_cartesian(ylim = c(-6,2), expand = FALSE)
#Export as 5 x 20 inch pdf

### STATS ###
TEDDY3to14vsTEDDY15to48 <- kruskal.test(Percent ~ Project, data=justTEDDY)
TEDDY3to14vsTEDDY15to48

TEDDYvsHMP <- subset(AllWGSpercent, Project!="NonWestern")
kruskal.test(Percent ~ Study, data=TEDDYvsHMP)

TEDDYvsNW <- subset(AllWGSpercent, Study!="HMP")
kruskal.test(Percent ~ Study, data=TEDDYvsNW)

#can also use wilcox.test