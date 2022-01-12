library(ggplot2)
library(RColorBrewer)
library(gridExtra)

colnames(fig02b)=c('source', 'ai', 'tavg')
display.brewer.pal(4,'Set2')
myColors <- brewer.pal(4,"Set2")

fig=ggplot(fig02b, aes(tavg, ai, color=source)) +
  geom_point(size=3) +
  theme_bw() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.075,0.80))+
  labs(y="Global Aridity Index", x=("Mean annual temperature (°C)")) +
  scale_colour_manual(name= 'Data source:', labels=c('aridec','WorldClim2'), values=c('#E78AC3','#8DA0CB')) +
  geom_hline(yintercept=c(0.05, 0.2, 0.5, 0.65), linetype= "dashed") +
  scale_y_continuous(breaks=c(0.05, 0.2, 0.5, 0.65, 1.0, 1.4))

ggsave("~/Documents/fig.pdf", device = "pdf")
