library("ggplot2"); packageVersion("ggplot2")
library(scales)

######### Supplementary Figure 2a: Top 5 Taxa  ##########
taxalist <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig2a.txt", header = TRUE, sep = "\t")

ggplot(taxalist, aes(Taxa, RelAbund, colour=Project)) +
  geom_point(aes(color=Project, shape=Project),size=8,position=position_jitterdodge(jitter.width=0.7,dodge.width = 0.8)) +
  theme_bw() +
  theme(
    text = element_text(size=40),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line    = element_line(color='black', size = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "black", size = 2),
    #plot.margin=margin(l=10,b=10,t=10,r=10,unit = "pt"),
    legend.position = "none",
    legend.title = element_blank()) +
  #scale_shape_manual(values=c(16, 2))+
  scale_fill_manual(values=c("red4","royalblue3")) +
  scale_colour_manual(values=c("red4","royalblue3")) +
  scale_y_continuous(trans='log10',limits=c(0.01,110), breaks=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100), labels = comma) +
  coord_cartesian(ylim = c(0.01,110), expand = FALSE)

#EXPORT FIGURE AS 4000 x 2000 tiff

######### Supplementary Figure 2b: Low p-values  ##########
taxalist <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig2b.txt", header = TRUE, sep = "\t")

ggplot(taxalist, aes(Taxa, RelAbund, colour=Project)) +
  geom_point(aes(color=Project, shape=Project),size=8,position=position_jitterdodge(jitter.width=0.7,dodge.width = 0.8)) +
  theme_bw() +
  theme(
    text = element_text(size=40),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line    = element_line(color='black', size = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "black", size = 2),
    #plot.margin=margin(l=10,b=10,t=10,r=10,unit = "pt"),
    legend.position = "none",
    legend.title = element_blank()) +
  scale_fill_manual(values=c("red4","royalblue3")) +
  scale_colour_manual(values=c("red4","royalblue3")) +
  scale_y_continuous(trans='log10',limits=c(0.01,110), breaks=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100), labels = comma) +
  coord_cartesian(ylim = c(0.01,110), expand = FALSE)

#EXPORT FIGURE AS 4000 x 2000 tiff