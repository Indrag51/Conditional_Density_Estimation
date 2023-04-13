library(ggplot2)

beta1<-theta.true
theta.true<- t(t(beta1)/sqrt(colSums(beta1 * beta1)))

beta<- t(t(beta)/sqrt(colSums(beta * beta)))

A = melt(theta.true)
names(A)<- c('coordinate','location','value')
B = melt(beta)
names(B)<- c('coordinate','location','value')

p11 = ggplot(A, aes(x = coordinate, y = location, fill = value)) +
  geom_tile() + scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"), limits=c(-1,1)) +
  coord_fixed() + labs(title = 'True_Beta') + theme_classic() + theme(plot.title = element_text(hjust=0.5))

p22 = ggplot(B, aes(x = coordinate, y = location, fill = -value)) +
  geom_tile() + scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"), limits=c(-1,1)) +
  coord_fixed() + labs(title = 'estimated_Beta') + theme_classic() + theme(plot.title = element_text(hjust=0.5))

gridExtra::grid.arrange(p11,p22,ncol=2)
