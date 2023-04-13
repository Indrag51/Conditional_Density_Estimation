load("C:\\Users\\indra\\Dropbox\\Debdeep_Neuro\\Paper_final\\code\\real_data_all_function.Rdata")
XI_mat = run1$XI_mat.out
BETA_mat = run1$BETA_mat.out

n <- nrow(X) #no. of individuals
K <- ncol(y) #no. of Location
p <- ncol(X) #no. of covariates

grid<- seq(0,1,length.out = 10001)
par(mfrow=c(2,2))
location.1<- c(1,10,20,30,40,50,60,70,80)
density_mat<- matrix(0,nr=length(grid),nc=length(location.1))
density_mat.true<- density_mat
density_mat.NNKCDE<- density_mat
ind.index = sort(sample(1:nrow(X_test),nrow(X_test),replace = F))
X_new=matrix(X_test[ind.index,],nc=p)
for(i in 1:length(location.1))
{
  density.array<- PostFun_loc(XI_mat, grid, X_new, BETA_mat,
                              location = location.1[i],test_knots = test_knots)
  for(j in 1:nrow(X_new)){
    density<- density.array[,,j]
    density_q2<- apply(density,1,function(x) {mean(x)})
    density_q2<- density_q2/sum(density_q2*(grid[2]-grid[1]))
    density_mat[,i] = density_mat[,i] + density_q2
    print(c(j,i))
  }
  density_mat[,i] = density_mat[,i]/nrow(X_new)
} 
par(mfrow=c(1,1))
dat1 = data.frame(grid,density_mat)
dat2 = data.frame(y_test[,location.1])
par(mfrow=c(2,2))
i=1
for(i in 1:length(location.1)){
  ggplot(dat1,aes(x=grid,y=density_mat[,i]),col="red") + geom_line(size=2,col="red") + xlab("\n grid") + ylab("Density \n") + 
    ggtitle(paste("Density at location ",location.1[i])) + geom_point(data = dat2, 
                                                                      aes(x=y_test[,location.1[i]],y=0),col="blue")
  i=i+1
}
par(mfrow=c(1,1))
