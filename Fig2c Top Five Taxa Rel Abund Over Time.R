library("ggplot2"); packageVersion("ggplot2")

taxa <- read.table(file="/Users/tommyauchtung/Desktop/Fig2c.txt", header = TRUE, sep = "\t")

ggplot(taxa) +
  theme_bw() +
  geom_smooth(method="loess", size=1, aes(x=ExactAgeinMonths, y=NinetyNine_Saccharomycescerevisiae1, colour = "Saccharomyces cerevisiae")) + 
  geom_smooth(method="loess", size=1, aes(x=ExactAgeinMonths, y=NinetyNine_Penicilliumpaneum2, colour = "Penicillium paneum")) +
  geom_smooth(method="loess", size=1, aes(x=ExactAgeinMonths, y=NinetyNine_Candidaalbicans3, colour = "Candida albicans")) +
  geom_smooth(method="loess", size=1, aes(x=ExactAgeinMonths, y=NinetyNine_Candidaparapsilosis4, colour = "Candida parapsilosis")) +
  geom_smooth(method="loess", size=1, aes(x=ExactAgeinMonths, y=NinetyNine_Candidazeylanoides5, colour = "Candida zeylanoides")) +
  labs(x = "Age (Months)", y = "Relative Abundance") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line    = element_line(color='black'),
        plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) +
  scale_colour_manual("", breaks = c("Saccharomyces cerevisiae", "Penicillium paneum", "Candida albicans", "Candida parapsilosis", "Candida zeylanoides"), values = c("yellow","forestgreen","saddlebrown","deeppink","dodgerblue2")) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(3,48)) +
  coord_cartesian(ylim = c(0, 0.4), xlim = c(0, 48), expand = FALSE) 

#Export image as 12 x 6 pdf 