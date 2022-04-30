library("ggplot2"); packageVersion("ggplot2")

#################### Figure 1a ###########################
percentfungiinWGS <- read.table(file="/Users/tommyauchtung/Desktop/Fig1a.txt", header = TRUE, sep = "\t")
#Percent Fungi Over Time
ggplot(percentfungiinWGS) +
  geom_smooth(aes(x=ExactAgeinMonths, y=LogPercent, colour = Country), method="loess", size = 4) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        #legend.title = element_blank(),
        #legend.position = "top",
        #unit(4, "cm"),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=20,b=15,t=15,r=0,unit = "pt")) +
  labs(x = "Age (Months)", y = "Percent Fungi") +
  scale_y_continuous(breaks=c(-5,-4,-3)) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 48), ylim = c(-5,-3), expand = FALSE) +
  scale_colour_manual(values=c("#002f6c", "#FF0000", "#FFCD00", "lawngreen"))

#Export 20 x 20 inch pdf
#Separately export image including country legend as large SVG image, then add/edit in Powerpoint

########## Supplementary Figure 4a: TOTAL FUNGI BY BREASTFEEDING STATUS #################
md <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig4a.txt", header = TRUE, sep = "\t")

ggplot(md) +
  geom_smooth(aes(x=ExactAgeinMonths, y=LogWGSpercent, colour = Breastfeeding), method="loess", size = 4) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(6, "cm"),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        text = element_text(size=30),
        legend.text=element_text(size=40),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=15,b=15,t=15,r=15,unit = "pt")) +
  labs(x = "Age (Months)", y = "Percent Fungi") +
  scale_y_continuous(breaks=c(-5,-4,-3)) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 48), ylim = c(-5,-3), expand = FALSE) +
  scale_colour_manual(values=c("mediumorchid4", "mediumseagreen"))

#Export 2000 x 2000 tiff

########## Supplementary Figure 4b,c,d: TOTAL FUNGI BY OTHER VARIABLES = PROBIOTICS, SEX, ASTHMA STATUS #################
md <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig4bcd.txt", header = TRUE, sep = "\t")
#Plot 'EverProbioticsYesNo' or 'Sex' or 'Asthma'
ggplot(md) +
  geom_smooth(aes(x=ExactAgeinMonths, y=LogWGSpercent, colour = EverProbioticsYesNo), method="loess", size = 4) +
  theme_bw() +
  theme(#axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        #axis.title.y = element_blank(),
        #legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(6, "cm"),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        text = element_text(size=30),
        legend.text=element_text(size=40),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=15,b=15,t=15,r=15,unit = "pt")) +
  labs(x = "Age (Months)", y = "Percent Fungi") +
  scale_y_continuous(breaks=c(-5,-4,-3)) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 48), ylim = c(-5,-3), expand = FALSE) +
  scale_colour_manual(values=c("mediumorchid4", "mediumseagreen"))

#Export 2000 x 2000 tiff
