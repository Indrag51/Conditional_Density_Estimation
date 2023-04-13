source("/Users/indrajit_ghosh/Dropbox/Debdeep_Neuro/Paper_final/code/Github/model_main/functions_main.R")
source("/Users/indrajit_ghosh/Dropbox/Debdeep_Neuro/Paper_final/code/Github/model_main/WAIC.R")
library(ddalpha)
#devtools::install_github("tpospisi/NNKCDE/r")
#library(NNKCDE)
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
K1=20
n2 = 0.2 * n1 ## test set
RNGkind(sample.kind = "Rejection")
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

test_knots = seq(0,1, length.out = 51)
Nplus = length(test_knots)

repetition<-50
TT<-vector("list",repetition)

for(tt in 1:repetition){
  set.seed(842*tt)
  X1<- matrix(runif(n1*p1,-5,5),nc=p1)
  X1<- apply(X1,MARGIN = 2,FUN = function(x)(x-mean(x)))
  max_norm<- max(sqrt(rowSums(X1*X1)))
  X2<- X1/max_norm
  
  beta<-theta.true
  theta.true<- t(t(beta)/sqrt(colSums(beta * beta)))
  W_0 <-(kronecker(diag(ncol(beta)),X2) %*% matrix(as.vector(theta.true),nc=1) + 1)/2
  W_1 <- 2*(W_0-0.5) #Centralize the probability
  y_vec<-numeric(nrow(X2) * ncol(theta.true))
  y1<- matrix(0,nr=n1,nc=K1)
  for (j in 1:K1) {
    for(i in 1:n1){
      # y_vec[(j-1)*nrow(X2) + i]<- rbeta(1,exp(2*W_0[(j-1)*nrow(X2) + i] + 1), 8*(2.25-j/ncol(theta.true)))
      y_vec[(j-1)*nrow(X2) + i]<- ifelse(runif(1)<=1/(1 + exp(-W_1[(j-1)*nrow(X2) + i])),rbeta(1,16/(1 + exp(W_0[(j-1)*nrow(X2) + i])),32*(j/ncol(theta.true))^0.5),
                                         rbeta(1,32*(sin(j/ncol(theta.true)))^{3/4},8*W_0[(j-1)*nrow(X2) + i]))
      }
  }
  
  y1<- matrix(y_vec,ncol = ncol(theta.true))
  X = X_train <- X2[-index.sample,]
  X_test <- X2[index.sample,]
  y = y_train <- y1[-index.sample,]
  y_test <- y1[index.sample,]
  
  beta_Sigma<- 1*Sigma_kron(tau=0.05,rho=0.75, K=ncol(y), p=ncol(X))
  
  Gamma1 = Matern_cov(test_knots, nu= 2.5, sigma = 1,theta = 0.292)
  
  CHOL1 = t(chol(Gamma1))
  xi = matrix(CHOL1%*%rnorm(Nplus^2),nc=Nplus,byrow = T)
  #xi = xi-mean(xi)
  beta<- matrix(as.vector(mvtnorm::rmvnorm(1,mean = rep(0, p1*K1))),nc=K1)
  beta<- t(t(beta)/sqrt(colSums(beta * beta)))
  tau = matrix(CHOL1%*%rnorm(Nplus^2),nc=Nplus,byrow = T)
  #tau = tau-mean(tau)
  h_loc = hBasis(seq(0,(K1-1),by=1)/(K1-1),test_knots)
  
  #######################
  #######################
  nmcmc=5000
  burnin=5000
  thining=5
  
  t0 = proc.time()
  run1<- get_updates(xi=xi,beta = beta, tau=tau,h_loc=h_loc, X=X,y=y,
                     test_knots = test_knots,
                     verbose = T, burnin = burnin,nmcmc = nmcmc,thining = thining)
  print(proc.time()-t0)
  
  
  print(paste("replication=",tt))
  
  XI_mat = run1$XI_mat.out
  TAU_mat = run1$TAU_mat.out
  BETA_mat = run1$BETA_mat.out
  
  n <- nrow(X) #no. of individuals
  K <- ncol(y) #no. of Location
  p <- ncol(X) #no. of covariates
  #######################
  #######################
  grid<- seq(0,1,length.out = 10001)
  par(mfrow=c(3,3))
  #location.1<-c(1,5,10,15,20)
  location.1<-sort(sample(1:K,K,replace = F))
  density_mat<- matrix(0,nr=length(grid),nc=length(location.1))
  density_mat.true<- density_mat
  #density_mat.NNKCDE<- density_mat
  #ind.index = c(9,11,12)
  ind.index = sort(sample(1:nrow(X_test),nrow(X_test),replace = F))
  X_new=matrix(X_test[ind.index,],nc=p)
  discrip1 = rep(0, length(ind.index))
  discrip2 = rep(0, length(location.1))
  
  for(i in 1:length(location.1))
  {
    theta.true.1<- theta.true[,location.1[i]]
    a.vec<- as.numeric((X_new%*%theta.true.1+1)/2)
    a1.vec<- (a.vec-0.5)*2
    density.array<- PostFun_loc(XI_mat, grid, X_new, BETA_mat, TAU_mat,
                                location = location.1[i],
                                test_knots = test_knots)
    for(j in 1:nrow(X_new)){
      a=a.vec[j]
      a1=a1.vec[j]
      #density.true<- dbeta(grid,exp(2*a+1), 8*(2.25-location.1[i]/ncol(theta.true)))
      density.true<- dbeta(grid,16/(1 + exp(a)),32*(location.1[i]/ncol(theta.true))^0.5)/(1+exp(-a1)) +
        dbeta(grid,32*(sin(location.1[i]/ncol(theta.true)))^{3/4},8*a)/(1+exp(a1))
      density<- density.array[,,j]
      density_q2<- apply(density,1,function(x) {quantile(x,0.50)})
      discrip1[j] = (sum((sqrt(density.true) - sqrt(density_q2))^2) * (range(grid)[2]-range(grid)[1])/length(grid)) ###Proxy hellinger distance
      if(j%%20==0) {print(paste("for individual =",j))}
    }
    discrip2[i] = mean(discrip1)
    if(i%%5==0) {print(paste("at location =",i))}
  }
  discrip3 = mean(discrip2)
  # par(mfrow=c(1,1))
  # plot(run1$LOG_Like.out,type='l', col='blue',lwd=2,ylab='Loglikelihood',xlab='Iteration')
  ####################################
  ####################################
  #save("run1", file = "bivariate_all_function_modified_3_implementation.Rdata")
  X_pred=X_test
  n=nrow(X_pred)
  Niter=nmcmc/thining
  Median_mat<- array(0,dim = c(n,K,Niter))
  u<- seq(0.05,0.95,by=0.05)
  #u<-c(0.15,0.30,0.45,0.6,0.75,.90)
  b<-rep(0,length(u))
  x.grid<- seq(0,1,length.out = 100001)
  Median_mat.loc = array(0,dim = c(n,K,length(u)))
  Median_mat.loc.th = array(0,dim = c(n,K,length(u)))
  
  xi=apply(XI_mat,MARGIN = c(1,2),mean)
  tau=apply(TAU_mat,MARGIN = c(1,2),mean)
  beta=apply(BETA_mat,MARGIN = c(1,2),mean)
  for(j in 1: length(u)){
    
    Median_mat.loc[,,j] = get_quantile(xi = xi,tau = tau,beta = beta,
                                       test_knots = test_knots,
                                       prob = u[j],X_pred = X_pred)
  }
  #Median_mat.loc ## quantiles vary row-wise and locations vary column-wise 
  
  
  for (k in 1:n) {
    for(i in 1:K){
      theta.true.1<- theta.true[,i]
      a<- as.numeric((X_pred[k,]%*%theta.true.1 + 1)/2)
      a1<- (a-0.5)*2
      # A = pbeta(x.grid,exp(2*a+1), 8*(2.25-i/ncol(theta.true)))
      # A = pbeta(x.grid,10*exp((i/K - 0.5)/2),4*i^{3/4})*1/(1+exp(a1)) +
      #   pbeta(x.grid,32*sin(a),8*(i/K+1)^{1/2})*1/(1+exp(-a1))
      A = pbeta(x.grid,16/(1 + exp(a)),32*(i/K)^0.5)/(1+exp(-a1)) +
        pbeta(grid,32*(sin(i/K))^{3/4},8*a)/(1+exp(a1))
      for(j in 1:length(u)){b[j]<- (x.grid[min(which(A>u[j]))])}
      Median_mat.loc.th[k,i,] = b
    }
    if(k%%10==0){print(k)}
  }
  #Median_mat.loc.th
  
  WAIC.train = waic.eval(XI_mat, y_train, X_train, BETA_mat, TAU_mat, test_knots)
  
  WAIC.test = waic.eval(XI_mat, y_test, X_test, BETA_mat, TAU_mat, test_knots)
  
  TT[[tt]]<- list(run1,discrip3, Median_mat.loc, Median_mat.loc.th, WAIC.test, WAIC.train)
}
save("TT", file = "Model_comparison.Rdata")
