library(ggplot2)

A = melt(theta.true)
names(A)<- c('coordinate','location','value')
B = melt(beta)
names(B)<- c('coordinate','location','value')

p1 = ggplot(A, aes(x = coordinate, y = location, fill = value)) +
  geom_tile() + scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"), limits=c(-1,1)) +
  coord_fixed() + labs(title = 'True_Beta') + theme_classic() + theme(plot.title = element_text(hjust=0.5))

p2 = ggplot(B, aes(x = coordinate, y = location, fill = value)) +
  geom_tile() + scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"), limits=c(-1,1)) +
  coord_fixed() + labs(title = 'estimated_Beta') + theme_classic() + theme(plot.title = element_text(hjust=0.5))

gridExtra::grid.arrange(p1,p2,ncol=2)
