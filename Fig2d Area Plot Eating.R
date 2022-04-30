library(ggplot2)

food <- read.table(file="/Users/tommyauchtung/Desktop/Fig2d.txt", header = TRUE, sep = "\t")

ggplot(food, aes(x = ExactAgeinMonths, y=..count.., fill = BreastFeedingFormulaFood)) + 
  geom_density(position="fill") +
  scale_fill_manual(values=c("chocolate2", "yellowgreen", "mediumpurple2", "darkgoldenrod1",  "deeppink3", "mediumseagreen")) +
  scale_y_continuous(labels=scales::percent) + 
  scale_x_continuous(breaks=c(0,6,12,18,24,30,36,42,48), limits = c(3,48)) +
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line    = element_line(color='black'),
        axis.title = element_text(size=28),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        #legend.position = c(0.7, 0.3),
        #legend.key.size = unit(2, "cm"),
        #legend.title = element_blank(),
        plot.margin=margin(l=20,b=15,t=15,r=20,unit = "pt")) +
  coord_cartesian(xlim = c(0, 48), expand = FALSE)

#Export image as 12 x 6 pdf 