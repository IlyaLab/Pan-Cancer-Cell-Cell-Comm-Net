
# plots
library(ggplot2)
library(dplyr)
library(forcats)

load('scaffold_median_cut_march11.rda')


scaffold$Cell <- factor(scaffold$Cell,
       levels=names(sort(table(scaffold$Cell),
             decreasing=TRUE)))

scaffold %>%
  ggplot( aes(x=Cell, fill=Type, col=Type) ) +
  geom_bar(position = 'dodge') +
  coord_flip() +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


ggplot(data=unique(scaffold[,c('Cell', 'Samples')]), aes(x=Samples)) +
  geom_bar(position = 'dodge') +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  ylab('Per Cell type')


ggplot(data=unique(scaffold[,c('Cell', 'Wt', 'Samples')]), aes(x=Wt, fill=Samples)) +
  geom_histogram()
  ylab('Per Cell type')





