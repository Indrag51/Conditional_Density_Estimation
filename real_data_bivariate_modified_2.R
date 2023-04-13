source("bivariate_all_function_modified_2.R")
library(ddalpha)
#devtools::install_github("tpospisi/NNKCDE/r")
library(NNKCDE)
library(xgboost)
####################
####################
set.seed(1616)
Covariates = read.table("Covariates.txt")
Response = read.table("Response.txt")
# Covariates = read.table("C:\\Users\\TZKNYG\\Documents\\Density_est\\Covariates.txt")
# Response = read.table("C:\\Users\\TZKNYG\\Documents\\Density_est\\Response.txt")

X<-lapply(Covariates[2:214,], function(x) as.numeric(as.character(x)))
X<- matrix(unlist(X),nr=213)
y<- matrix(unlist(Response),nr=213)
#y<- y[,sample(1:83,size = 50,replace = F)] 
max_norm<- max(sqrt(rowSums(X*X)))
X2<- X/max_norm
y1=y
n1 = nrow(X)
n2 = 0.2 * n1 ## test set
RNGkind(sample.kind = "Rejection")
set.seed(12)
index.sample = sample(1:n1,size=n2,replace = F)
X = X_train <- X2[-index.sample,]
X_test <- X2[index.sample,]
y = y_train <- y1[-index.sample,]
y_test <- y1[index.sample,]

test_knots = seq(0,1, length.out = 51)
Nplus = length(test_knots)
K1=ncol(y)
beta_Sigma<- 4*Sigma_kron(tau=0.25,rho=0.25, K=ncol(y), p=ncol(X))
#test_knots2 = seq(0,1,length.out = 101)
Gamma1 = kronecker(CovMat.2(test_knots, nu= 2.5, sigma = 2,theta = .25),
                   CovMat.2(test_knots, nu= 2.5, sigma = 1,theta = .5))
# Gamma1 = kronecker(exp(-500*outer(1:Nplus,1:Nplus,FUN = function(x,y) {return(((x-y)/Nplus)^2)})),
#            CovMat.2(test_knots, nu= 2.5, sigma = 2,theta = .15))
CHOL1 = t(chol(Gamma1))
xi = matrix(CHOL1%*%rnorm(Nplus^2),nc=Nplus,byrow = T)
xi = xi-mean(xi)
beta<- matrix(as.vector(mvtnorm::rmvnorm(1,mean = rep(0, ncol(X)*ncol(y)))),nc=ncol(y))
beta<- t(t(beta)/sqrt(colSums(beta * beta)))
tau = matrix(CHOL1%*%rnorm(Nplus^2),nc=Nplus,byrow = T)
tau = tau-mean(tau)
h_loc = hBasis(seq(0,(K1-1),by=1)/(K1-1),test_knots)
#######################
#######################
nmcmc=5000
burnin=15000
thining=10
run1<- get_updates(xi=xi,beta = beta, tau=tau,h_loc=h_loc, X=X,y=y,
                   test_knots = test_knots,
                   verbose = T, burnin = burnin,nmcmc = nmcmc,thining = thining)

XI_mat = run1$XI_mat.out
TAU_mat = run1$TAU_mat.out
BETA_mat = run1$BETA_mat.out

n <- nrow(X) #no. of individuals
K <- ncol(y) #no. of Location
p <- ncol(X) #no. of covariates
#######################
#save("run1", file = "parameters_beta_grid_1001.RData")
#######################
grid<- seq(0,1,length.out = 10001)
par(mfrow=c(3,3))
#location.1<-sort(sample(1:K,K,replace = F))
location.1<- c(1,10,20,30,40,50,61,72,83)
density_mat<- matrix(0,nr=length(grid),nc=length(location.1))
density_mat.true<- density_mat
density_mat.NNKCDE<- density_mat
ind.index = sort(sample(1:nrow(X_test),3,replace = F))
X_new=matrix(X_test[ind.index,],nc=p)
for(i in 1:length(location.1))
{
  # theta.true.1<- theta.true[,location.1[i]]
  # a.vec<- as.numeric((X_new%*%theta.true.1+1)/2)
  # a1.vec<- (a.vec-0.5)*2
  density.array<- PostFun_loc(XI_mat, grid, X_new, BETA_mat, TAU_mat,
                              location = location.1[i],
                              test_knots = test_knots)
  for(j in 1:nrow(X_new)){
    #a=a.vec[j]
    #a1=a1.vec[j]
    # density.true<- dbeta(grid,10*exp((location.1[i]/ncol(theta.true) - 0.5)/2),30)*1/(1+exp(a1)) +
    #   dbeta(grid,32*sin(a),8*(location.1[i]/ncol(theta.true)+1)^{1/2})*1/(1+exp(-a1))
    #density.true<- dbeta(grid,8*a,40*(location.1[i]/ncol(theta.true))) * 0.6 + 0.4*dbeta(grid,10*exp(a-0.5),4)
    # dbeta(grid,3*location.1[i]^{3/4},8*a*(2-location.1[i]/K))*0.3)
    density<- density.array[,,j]
    density_q2<- apply(density,1,function(x) {quantile(x,0.50)})
    density_q2<- density_q2/sum(density_q2*(grid[2]-grid[1]))
    density_q1<- apply(density,1,function(x) {quantile(x,0.025)})
    density_q3<- apply(density,1,function(x) {quantile(x,0.975)})
    max_lim<- min(max(density_q2))
    min_lim<- min(density_q2)
    plot(grid, density_q2,ylab = "",ylim= c(0, 10), type = "l",col="red",lwd=2,lty=1,
         main = paste("Overlay plot at loc=",location.1[i],", individual =",ind.index[j]))
    #lines(grid, density_q2,type = "l",col="blue",lwd=2,lty=1)
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
save("run1", file = "real_data_modified_2.Rdata")
