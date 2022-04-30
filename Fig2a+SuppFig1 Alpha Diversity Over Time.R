library("ggplot2"); packageVersion("ggplot2")

alphadiversity <- read.table(file="/Users/tommyauchtung/Desktop/Fig2a+SuppFig1.txt", header = TRUE, sep = "\t")

#################### Figure 1a ###########################
#16S + ITS Alpha Diversity Over Time (SHANNON)
ggplot(alphadiversity) +
  geom_smooth(aes(x=ExactAgeinMonths, y=Shan16SR3K), method="loess", size = 1, color = "red") +
  geom_smooth(aes(x=ExactAgeinMonths, y=ShanFungiR3K), method="loess", size = 1, color = "blue") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) +
  labs(x = "Age (Months)", y = "Shannon Diversity") +
  coord_cartesian(xlim = c(0, 48), ylim = c(0, 3.25), expand = FALSE) +
  geom_point(size=0.1, (aes(x=ExactAgeinMonths, y=Shan16SR3K)), color = "red") +
  geom_point(size=0.1, (aes(x=ExactAgeinMonths, y=ShanFungiR3K)), color = "blue") +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48)) 

#Export as 12 x 6 pdf

#################### Supplementary Figure 1 ###########################
#16S + ITS Alpha Diversity Over Time (OBSERVED)
ggplot(alphadiversity) +
  geom_smooth(aes(x=ExactAgeinMonths, y=Obs16SR3K), method="loess", size = 2, color = "red") +
  geom_smooth(aes(x=ExactAgeinMonths, y=ObsFungiR3K), method="loess", size = 2, color = "blue") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=margin(l=15,b=15,t=15,r=25,unit = "pt")) +
  labs(x = "Age (Months)", y = "Observed Diversity") +
  geom_point(size=0.5, (aes(x=ExactAgeinMonths, y=Obs16SR3K)), color = "red") +
  geom_point(size=0.5, (aes(x=ExactAgeinMonths, y=ObsFungiR3K)), color = "blue") +
  coord_cartesian(xlim = c(0, 48), ylim = c(0, 100), expand = FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48))

#Export as 2000 x 1000 tiff

### STATS ###
TEDDY3to15vsTEDDY16to48 <- subset(alphadiversity, ExactAgeinMonths > 2 & ExactAgeinMonths < 49)
#Fungi
wilcox.test(ShanFungiR3K ~ RoughAge, data=TEDDY3to15vsTEDDY16to48)
wilcox.test(ObsFungiR3K ~ RoughAge, data=TEDDY3to15vsTEDDY16to48)
#Bacteria
wilcox.test(Shan16SR3K ~ RoughAge, data=TEDDY3to15vsTEDDY16to48)
wilcox.test(Obs16SR3K ~ RoughAge, data=TEDDY3to15vsTEDDY16to48)

