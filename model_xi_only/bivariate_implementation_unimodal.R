source("bivariate_all_function.R")
library(ddalpha)
#devtools::install_github("tpospisi/NNKCDE/r")
library(NNKCDE)
library(xgboost)
# devtools::install_github("rizbicki/FlexCoDE")
# library(FlexCoDE)
#devtools::install_github("tpospisi/cdetools/r")
#library(cdetools)
##################################################
##### Data Generation & Initialization ###########
##################################################
n1=200
p1=5
K1=30
n2 = 0.2 * n1 ## test set
set.seed(403)
index.sample = sample(1:n1,size=n2,replace = F)
theta_func<- function(j,x){
  if(j==1)
    return(x^2+1)
  if(j==2)
    return((1-x)^2)
  if(j==3)
    return(4*x*(1-x))
  if(j==4)
    return(-1 + 2*(x-0.75)^2)
  if(j==5)
    return(1.5 + 4*sin(x-0.5))
  if(j==6)
    return(cos(x+2) + 0.25 * (x-0.5)^4 )
  if(j==7)
    return(log(x+1)- 4*x^2)
  else
    return("wrong input")
}
theta.true<-matrix(0,nr=p1,nc=K1)
for(j in 1:p1){
  theta.true[j,]<- (matrix(theta_func(j, seq(0,1,length.out = K1)),nr=1))
}

set.seed(842)
X1<- matrix(runif(n1*p1,-5,5),nc=p1)
X1<- apply(X1,MARGIN = 2,FUN = function(x)(x-mean(x)))
max_norm<- max(sqrt(rowSums(X1*X1)))
X2<- X1/max_norm

#X_new<- kronecker(diag(K1),X)
# theta.true<- matrix(rnorm(p*K),nr=p)
beta<-theta.true
theta.true<- t(t(beta)/sqrt(colSums(beta * beta)))
W_0 <-(kronecker(diag(ncol(beta)),X2) %*% matrix(as.vector(theta.true),nc=1) + 1)/2
W_1 <- 2*(W_0-0.5) #Centralize the probability
y_vec<-numeric(nrow(X2) * ncol(theta.true))
for (j in 1:ncol(theta.true)) {
  for(i in 1:nrow(X2)){
    y_vec[(j-1)*nrow(X2) + i]<-ifelse(runif(1)<=1,rbeta(1,16*W_0[(j-1)*nrow(X2) + i],32),
                                      rbeta(1,16*W_0[(j-1)*nrow(X2) + i],16))
    
    # y_vec[(j-1)*nrow(X2) + i]<-ifelse(runif(1)<=1/(1+exp(W_1[(j-1)*nrow(X2) + i])),
    #                                   rbeta(1,10*exp((j/ncol(theta.true)-0.5)/2),30),
    #                                   rbeta(1,32*sin(W_0[(j-1)*nrow(X2) + i]), 8*(j/ncol(theta.true) +1)^{1/2}))
    
    # y_vec[(j-1)*nrow(X) + i]<- ifelse(runif(1)<=0.5,rbeta(1,16*W_0[(j-1)*nrow(X) + i],16*j/ncol(theta.true)),
    #                                   rbeta(1,16,16*W_0[(j-1)*nrow(X) + i]))
    # y_vec[(j-1)*nrow(X) + i]<- ifelse(runif(1)<=0.5,rbeta(1,16*W_0[(j-1)*nrow(X) + i],32*j/ncol(theta.true)),
    #                                   rbeta(1,32,16*W_0[(j-1)*nrow(X) + i]))
    # y_vec[(j-1)*nrow(X) + i]<- ifelse(runif(1)<=0.5 ,rbeta(1,4*W_0[(j-1)*nrow(X) + i],5*abs(sin(32*j/ncol(theta.true)))),
    #                                   rbeta(1,(3*j/ncol(theta.true))^{3/4},8*W_0[(j-1)*nrow(X) + i]*(2-j/ncol(theta.true))))
    
  }
}
y1<- matrix(y_vec,ncol = ncol(theta.true))
test_knots = seq(0,1, length.out = 51)
Nplus = length(test_knots)
X = X_train <- X2[-index.sample,]
X_test <- X2[index.sample,]
y = y_train <- y1[-index.sample,]
y_test <- y1[index.sample,]

beta_Sigma<- 4*Sigma_kron(tau=0.25,rho=0.25, K=ncol(y), p=ncol(X))
#test_knots2 = seq(0,1,length.out = 101)
Gamma1 = kronecker(CovMat.2(test_knots, nu= 4.5, sigma = 0.5,theta = .10),
                   CovMat.2(test_knots, nu= 2.5, sigma = 2,theta = .75))
# Gamma1 = kronecker(exp(-500*outer(1:Nplus,1:Nplus,FUN = function(x,y) {return(((x-y)/Nplus)^2)})),
#            CovMat.2(test_knots, nu= 2.5, sigma = 2,theta = .15))
CHOL1 = t(chol(Gamma1))
xi = matrix(CHOL1%*%rnorm(Nplus^2),nc=Nplus,byrow = T)
xi = xi-mean(xi)
beta<- matrix(as.vector(mvtnorm::rmvnorm(1,mean = rep(0, ncol(X)*ncol(y)),sigma = beta_Sigma)),nc=ncol(y))
beta<- t(t(beta)/sqrt(colSums(beta * beta)))
#######################
#######################
nmcmc=2500
burnin=2500
thining=5
run1<- get_updates(xi=xi,beta = beta,X=X,y=y,
                   test_knots = test_knots,
                   verbose = T, burnin = burnin,nmcmc = nmcmc,thining = thining)

XI_mat = run1$XI_mat.out
BETA_mat = run1$BETA_mat.out

n <- nrow(X) #no. of individuals
K <- ncol(y) #no. of Location
p <- ncol(X) #no. of covariates
#######################
#######################
set.seed(738)
grid<- seq(0,1,length.out = 10001)
par(mfrow=c(3,3))
#location.1<-sort(sample(1:K,3,replace = F))
location.1<- c(1,9,15)
density_mat<- matrix(0,nr=length(grid),nc=length(location.1))
density_mat.true<- density_mat
density_mat.NNKCDE<- density_mat
ind.index = sample(1:nrow(X_test),3,replace = F)
X_new=matrix(X_test[ind.index,],nc=p)
for(i in 1:length(location.1))
{
  theta.true.1<- theta.true[,location.1[i]]
  a.vec<- as.numeric((X_new%*%theta.true.1+1)/2)
  a1.vec<- (a.vec-0.5)*2
  density.array<- PostFun_loc(XI_mat, grid, X_new, BETA_mat,location = location.1[i],
                              test_knots = test_knots, burnin = burnin,nmcmc = nmcmc,thining = thining)
  for(j in 1:nrow(X_new)){
    a=a.vec[j]
    a1=a1.vec[j]
    # density.true<- dbeta(grid,10*exp((location.1[i]/ncol(theta.true) - 0.5)/2),30)*1/(1+exp(a1)) +
    #   dbeta(grid,32*sin(a),8*(location.1[i]/ncol(theta.true)+1)^{1/2})*1/(1+exp(-a1))
    density.true<- (dbeta(grid,16*a,32) * 4 +
                      dbeta(grid,16*a,16) * 0)/4
    density<- density.array[,,j]
    density_q2<- apply(density,1,function(x) {quantile(x,0.50)})
    density_q2<- density_q2/sum(density_q2*(grid[2]-grid[1]))
    density_q1<- apply(density,1,function(x) {quantile(x,0.025)})
    density_q3<- apply(density,1,function(x) {quantile(x,0.975)})
    max_lim<- min(max(density.true,density_q2),10^3)
    min_lim<- min(density.true,density_q2)
    plot(grid, density.true,ylab = "",ylim= c(0, 7.3), type = "l",col="red",lwd=2,lty=1,
         main = paste("Overlay plot at loc=",location.1[i],", individual =",ind.index[j]))
    lines(grid, density_q2,type = "l",col="blue",lwd=2,lty=1)
    lines(grid, density_q1,type = "l",col="purple",lwd=2,lty=4)
    lines(grid, density_q3,type = "l",col="purple",lwd=2,lty=4)
    print(i)
    #print(sum((sqrt(density.true) - sqrt(density_q2))^2) * (range(grid)[2]-range(grid)[1])/length(grid)) ###Proxy hellinger distance
    # density_mat[,i]<- density_q2
    # density_mat.true[,i]<-density.true
  }
  # M_pred = matrix(c(matrix(X[d.q,],nc=p),location.1[i]),nr=1)
  # fit_pred = fit.NNKCDE$predict(M_pred,grid,h=0.1,k=45)
  # fit_pred = fit_pred/sum(fit_pred*(grid[2]-grid[1]))
  # density_mat.NNKCDE[,i] = matrix(fit_pred,nc=1)
  # lines(grid, density_mat.NNKCDE[,i],type = "l",col="green",lwd=2,lty=1)
} 
par(mfrow=c(1,1))
plot(run1$LOG_Like.out,type='l', col='blue',lwd=2,ylab='Loglikelihood',xlab='Iteration')
####################################