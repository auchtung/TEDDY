library("ggplot2"); packageVersion("ggplot2")

#################### Supplementary Figure 5 ###########################
regionaldiffs <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig5.txt", header = TRUE, sep = "\t")

# Supplementary Fig 5a: Clavispora lusitaniae over time (higher in early Florida/Georgia)
ggplot(regionaldiffs, aes(x=ExactAgeinMonths, y=Clavispora_lusitaniae_6, colour = Clinical.Center)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_color_brewer(palette="Set1") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) + 
  coord_cartesian(ylim = c(0, 0.1501,0.05), expand = FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48)) +
  labs(x = "Age (Years)", y = "Relative Abundance of Clavispora lusitaniae")

# Supplementary Fig 5b: Candida_zeylanoides over time (higher in Florida/Georgia)
ggplot(regionaldiffs, aes(x=ExactAgeinMonths, y=Candida_zeylanoides_5, colour = Clinical.Center)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_color_brewer(palette="Set1") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) + 
  coord_cartesian(ylim = c(0, 0.1501,0.05), expand = FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48)) +
  labs(x = "Age (Years)", y = "Relative Abundance of Candida_zeylanoides")

# Supplementary Fig 5c: Geotrichum candidum over time (higher in Finland)
ggplot(regionaldiffs, aes(x=ExactAgeinMonths, y=Galactomyces_candidum_7, colour = Clinical.Center)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_color_brewer(palette="Set1") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) + 
  coord_cartesian(ylim = c(0, 0.1501,0.05), expand = FALSE) +
  labs(x = "Age (Months)", y = "Relative Abundance of Geotrichum candidum") +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48))

#################### Supplementary Figure 6 ###########################
JustT1D <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig6.txt", header = TRUE, sep = "\t")

# Supplementary Fig 6a: Candida_albicans over time
ggplot(JustT1D, aes(x=ExactAgeinMonths, y=Candida_albicans_3, colour = T1D.OutcomeYesNoCorrected)) +
  geom_smooth(method="loess", size=8) + 
  theme_bw() +
  scale_colour_manual(values=c("seagreen3", "blueviolet")) + 
  theme(axis.text=element_text(size=26), 
        axis.title=element_text(size=26),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        #legend.position = "right",
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) +
  #legend.key.size = unit(1, "in"))  +
  labs(x = "Age (Months)", y = "Relative Abundance of Candida_albicans") +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48)) +
  coord_cartesian(ylim = c(0, 0.2, 0.05), expand = FALSE)

# Supplementary Fig 6b: Candida over time
ggplot(JustT1D, aes(x=ExactAgeinMonths, y=Candida, colour = T1D.OutcomeYesNoCorrected)) +
  geom_smooth(method="loess", size=8) + 
  theme_bw() +
  scale_colour_manual(values=c("seagreen3", "blueviolet")) + 
  theme(axis.text=element_text(size=26), 
        axis.title=element_text(size=26),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        #legend.position = "right",
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) +
  #legend.key.size = unit(1, "in"))  
  coord_cartesian(ylim = c(0, 0.4, 0.1), expand = FALSE) +
  labs(x = "Age (Months)", y = "Relative Abundance of Candida") +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48))

#################### Supplementary Figure 7 ###########################
JustIAPreCeliac <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig7.txt", header = TRUE, sep = "\t")

# Supplementary Fig 7a: Candida sake over time (higher in Celiac)
ggplot(JustIAPreCeliac, aes(x=ExactAgeinMonths, y=Candida_sake_18, color=UndiagnosedCeliac)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_colour_manual(values=c("dodgerblue2","orange")) + 
  labs(x = "Age (Years)", y = "Relative Abundance of Candida sake") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        #legend.position = "right",
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) +
  #legend.key.size = unit(1, "in")) + 
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(3,48)) +
  scale_y_continuous(breaks = seq(0,0.05, by=0.01)) +
  coord_cartesian(ylim = c(0, 0.05, 0.01),xlim = c(0, 48), expand = FALSE) 

# Supplementary Fig 7b: Candida albicans over time (NOT higher in Celiac)
ggplot(JustIAPreCeliac, aes(x=ExactAgeinMonths, y=Candida_albicans_3, color=UndiagnosedCeliac)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_colour_manual(values=c("dodgerblue2","orange")) + 
  labs(x = "Age (Years)", y = "Relative Abundance of Candida albicans") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        #legend.position = "right",
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) +
  #legend.key.size = unit(1, "in")) + 
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(3,48)) +
  scale_y_continuous(breaks = seq(0,0.2, by=0.05)) +
  coord_cartesian(ylim = c(0, 0.2, 0.05),xlim = c(0, 48), expand = FALSE) 

#################### Supplementary Figure 8 ###########################
Ppaneum <- read.table(file="/Users/tommyauchtung/Desktop/SuppFig8.txt", header = TRUE, sep = "\t")

# Supplementary Figure 8a: P.paneum over time (By Country)
ggplot(Ppaneum, aes(x=ExactAgeinMonths, y=Penicillium_paneum_2, colour = Clinical.Country)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position = "none",
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) + 
  coord_cartesian(ylim = c(0, 0.401,0.1), expand = FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48)) +
  scale_colour_manual(values=c("#336699", "#993399","#FF6600", "lawngreen")) + 
  labs(x = "Age (Years)", y = "Relative Abundance of P.paneum")

# Supplementary Figure 8b: P.paneum over time (By Breastfeeding)
ggplot(Ppaneum, aes(x=ExactAgeinMonths, y=Penicillium_paneum_2, colour = BreastfeedingYesNo)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_colour_manual(values=c("violetred1", "skyblue3")) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) + 
  coord_cartesian(ylim = c(0, 0.401,0.1), expand = FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48)) +
  labs(x = "Age (Years)", y = "Relative Abundance of P.paneum")

# Supplementary Figure 8c: P.paneum over time (By Formula)
ggplot(Ppaneum, aes(x=ExactAgeinMonths, y=Penicillium_paneum_2, colour = OnFormula)) +
  geom_smooth(method="loess", size=8) +
  theme_bw() +
  scale_color_brewer(palette="Dark2") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line    = element_line(color='black',size = 1, linetype = "solid"),
        axis.ticks.length = unit(.3, "cm"),
        axis.ticks = element_line(size = 1)) + 
  coord_cartesian(ylim = c(0, 0.401,0.1), expand = FALSE) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(0,48)) +
  labs(x = "Age (Years)", y = "Relative Abundance of P.paneum")

#Export images as 2000 x 2000 tiffs 


### SUPPLEMENTARY FIGURE 7b STATS ###
#To calculate C. albicans stats, need additional information
JustT1D <- read.table(file="/Users/tommyauchtung/Desktop/taxa+metadatatogether_JustT1D.txt", header = TRUE, sep = "\t")

### Subset by age for analysis of all time windows ###
JustT1D3to6 <- subset(JustT1D, ExactAgeinMonths >= 3 & ExactAgeinMonths < 7)
JustT1D3to6 <- plyr::ddply(JustT1D3to6, 'SubjectMaskId', head, n=1)
JustT1D7to10 <- subset(JustT1D, ExactAgeinMonths >= 7 & ExactAgeinMonths < 11)
JustT1D7to10 <- plyr::ddply(JustT1D7to10, 'SubjectMaskId', head, n=1)
JustT1D11to14 <- subset(JustT1D, ExactAgeinMonths >= 11 & ExactAgeinMonths < 15)
JustT1D11to14 <- plyr::ddply(JustT1D11to14, 'SubjectMaskId', head, n=1)
JustT1D15to18 <- subset(JustT1D, ExactAgeinMonths >= 15 & ExactAgeinMonths < 19)
JustT1D15to18 <- plyr::ddply(JustT1D15to18, 'SubjectMaskId', head, n=1)
JustT1D19to22 <- subset(JustT1D, ExactAgeinMonths >= 19 & ExactAgeinMonths < 23)
JustT1D19to22 <- plyr::ddply(JustT1D19to22, 'SubjectMaskId', head, n=1)
JustT1D23to26 <- subset(JustT1D, ExactAgeinMonths >= 23 & ExactAgeinMonths < 27)
JustT1D23to26 <- plyr::ddply(JustT1D23to26, 'SubjectMaskId', head, n=1)
JustT1D27to30 <- subset(JustT1D, ExactAgeinMonths >= 27 & ExactAgeinMonths < 31)
JustT1D27to30 <- plyr::ddply(JustT1D27to30, 'SubjectMaskId', head, n=1)
JustT1D31to34 <- subset(JustT1D, ExactAgeinMonths >= 31 & ExactAgeinMonths < 35)
JustT1D31to34 <- plyr::ddply(JustT1D31to34, 'SubjectMaskId', head, n=1)
JustT1D35to38 <- subset(JustT1D, ExactAgeinMonths >= 35 & ExactAgeinMonths < 39)
JustT1D35to38 <- plyr::ddply(JustT1D35to38, 'SubjectMaskId', head, n=1)
JustT1D39to42 <- subset(JustT1D, ExactAgeinMonths >= 39 & ExactAgeinMonths < 43)
JustT1D39to42 <- plyr::ddply(JustT1D39to42, 'SubjectMaskId', head, n=1)

wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D3to6)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D3to6)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D7to10)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D7to10)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D11to14)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D11to14)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D15to18)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D15to18)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D19to22)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D19to22)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D23to26)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D23to26)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D27to30)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D27to30)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D31to34)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D31to34)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D35to38)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D35to38)
wilcox.test(Candida ~ T1D.OutcomeYesNoCorrected, data=JustT1D39to42)
wilcox.test(Candida_albicans_3 ~ T1D.OutcomeYesNoCorrected, data=JustT1D39to42)