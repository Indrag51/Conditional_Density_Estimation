#library(truncnorm)
library(mvtnorm)
#library(spate)
#library(PrevMap)
library(fields)
library(ggplot2)
library(rSPDE)#Separable structure on Dispersion of beta
library(Rcpp)
#library(mcmcplots)
library(reshape2)
sourceCpp("bivariate_function_modified_2.cpp")

# get_sample = function(xi, tau, zeta, omega, X_new, beta, test_knots){
#   knot_length = length(test_knots)
#   intsize = 1/(knot_length-1)
#   p=nrow(beta)
#   K=ncol(beta)
#   n=nrow(X_new)/K
#   location_grid<- seq(0,K-1,by=1)/(K-1)
#   h_loc<-hBasis(location_grid,test_knots2)
#   g4<- h_loc%*%omega
#   beta_tilde <- t(t(beta)/sqrt(colSums(beta * beta)))
#   beta_tilde_aug<- matrix(as.vector(beta_tilde),nc=1)
#   Z<- (X_new%*%beta_tilde_aug + 1)/2
#   H_Z <- PhiBasis(Z,test_knots)
#   g2<- (zeta[1] + H_Z%*%zeta[2:length(zeta)])
#   delta<- matrix(xi,nc=1)%*%matrix(g2,nr=1) + matrix(tau,nc=1)%*%matrix(rep(g4,each=n),nr=1)
#   delta.max = max(delta)
#   A2 = delta[2:nrow(delta),]
#   A1 = delta[1:(nrow(delta)-1),]
#   
#   A3 = exp(delta.max)*(exp(A2-delta.max)-exp(A1-delta.max))/((A2-A1)/intsize)
#   
#   A4 = apply(A3,MARGIN = 2, function(x){cumsum(x)})
#   
#   constant = A4[nrow(A4),]
#   
#   u = runif(n*K)
#   
#   s=rep(0,n*K)
#   for(i in 1:(n*K)){
#     alpha = u[i]
#     const = constant[i]
#     j_0 = max(which(A4[,i]<alpha * const),0)+1
#     if(j_0>1)
#       s[i] = test_knots[j_0+1] + (log((alpha * const - A4[(j_0-1),i]) * 
#                                         (((A2-A1)/intsize)[j_0,i]) + exp(delta[j_0,i])) - 
#                                     delta[(j_0+1),i])/(((A2-A1)/intsize)[j_0,i])
#     else
#       s[i] = test_knots[j_0+1] + (log((alpha * const) * 
#                                         (((A2-A1)/intsize)[j_0,i]) + exp(delta[j_0,i])) - 
#                                     delta[(j_0+1),i])/(((A2-A1)/intsize)[j_0,i])
#   }
#   return(s)
# }
# 
# get_plot = function(y_grid = seq(0.0001,1,length.out = 1000), loc_test = 2, xi.truth,
#                     tau.truth,omega.truth,zeta.truth, theta.true,X_med, test_knots){
#   if(loc_test>K1) stop("Location is out of bounds!")
#   knot_length = length(test_knots)
#   intsize = 1/(knot_length-1)
#   
#   g4 = hBasis((loc_test-1)/(K1-1),test_knots2)%*%omega.truth
#   
#   if(missing(X_med))
#     X_med=matrix(apply(X, 2, function(x) {median(x)}),nc=p1)
#   
#   Z = (X_med%*%theta.true[,loc_test]+1)/2
#   g2 = zeta.truth[1] + PhiBasis(Z,test_knots)%*%zeta.truth[2:length(zeta.truth)]
#   delta<- matrix(xi.truth,nc=1)%*%matrix(g2,nr=1) + matrix(tau.truth,nc=1)%*%matrix(rep(g4,each=nrow(X_med)),nr=1)
#   delta.max = max(delta)
#   A2 = delta[2:nrow(delta),]
#   A1 = delta[1:(nrow(delta)-1),]
#   
#   A3 = exp(delta.max)*(exp(A2-delta.max)-exp(A1-delta.max))/((A2-A1)/intsize)
#   if(nrow(X_med)>1){
#     A4 = apply(A3,MARGIN = 2, function(x){sum(x)})
#     constant = sum(A4)
#   } else {
#     constant = sum(A3)
#   }
#   plot(y_grid, exp(hBasis(y_grid,test_knots)%*%delta)/constant,ylab="f(y|X,s)",xlab="grid")  
# }
# 
Sigma_kron<- function(tau = 1, rho = 0.5, K, p){
  H<- exp(-tau * outer(1:K,1:K,FUN = function(x,y) {return((x-y)^2)}))
  T_rho<- (1-rho) * diag(rep(1,p)) + rho * matrix(1,nc=p,nr=p)
  return(kronecker(H,T_rho))
}

# Matern kernel with smoothness parameter 3/2 or 5/2
mk = function(d, sigma=1, theta=1, flag=1)
{
  if(flag == 0){
    return((sigma^2)*(1+(sqrt(3)*abs(d)/theta))*exp(-sqrt(3)*abs(d)/theta))} # matern 3/2
  else{
    return((sigma^2)*(1+sqrt(5)*abs(d)/theta+5*(d)^2/(3*theta^2))*exp(-sqrt(5)*abs(d)/theta))} # matern 5/2
}
# first derivative of Matern kernel with respect to the first variable
mkp = function(d, sigma=1, theta=1, flag=1){
  if(flag == 0){
    return(-sigma^2*((3/theta^2)*abs(d)*exp(-sqrt(3)/theta*abs(d))))} # matern 3/2
  else{
    return(-(sigma^2)*(5/(3*theta^2)*abs(d) + 5*sqrt(5)*(d)^2/(3*theta^3))*exp(-sqrt(5)*abs(d)/theta))} # matern 5/2
}
# second derivative of Matern kernel with respect to the first and second variables
mkpp = function(d, sigma=1, theta=1, flag=1){
  if(flag == 0){
    return(sigma^2*((3/theta^2)*exp(-sqrt(3)/theta*abs(d))*(1-sqrt(3)/theta*abs(d))))} # matern 3/2
  else{
    return((sigma^2)*(5/(3*theta^2)+5*sqrt(5)*abs(d)/(3*theta^3)-25*(d)^2/(3*theta^4))*exp(-sqrt(5)*abs(d)/theta))} # matern 5/2
}

#####Covariance matrices for smoothness 3/2 and 5/2
CovMat.1 = function(knots, sigma = 1, theta = 1, flag = 1)
{
  Nknots = length(knots)
  CovMat = matrix(0, Nknots, Nknots)
  d<-as.matrix(dist(knots, method =  "euclidean", diag = T, upper = T))
  CovMat<- mk(d, sigma, theta, flag)
  return(CovMat)
}

#### This is for general smoothness ####
Matern_cov = function(test_knots, nu=2.5, sigma=1, theta=0.5){
  Nplus = length(test_knots)
  A = data.frame(x=rep(test_knots,each=Nplus),y=rep(test_knots,Nplus))
  B = (as.matrix(dist(A, method = "euclidean", diag = T, upper = T)))
  C = Matern(B,alpha = sqrt(2*nu)/theta, nu=nu, phi = sigma^2)
  return(C)
}

CovMat.2 = function(knots, sigma = 1, theta = 1, nu = 1.5)
{
  Nknots = length(knots)
  CovMat = matrix(0, Nknots, Nknots)
  d<-as.matrix(dist(knots, method =  "euclidean", diag = T, upper = T))
  CovMat<- Matern(d,alpha = sqrt(2*nu)/theta, nu=nu, phi = sigma^2)
  return(CovMat)
}

CovMat_monotone = function(knots, sigma = 1, theta = 1, flag = 1)
{
  Nknots = length(knots)
  CovMat = matrix(0, Nknots+1, Nknots+1)
  d<-as.matrix(dist(knots, method =  "euclidean", diag = T, upper = T))
  CovMat[1,1]<- mk(d[1,1], sigma, theta, flag)
  CovMat[1,2:(Nknots+1)]<- sapply(d[1,1:Nknots], FUN = function(x){mkp(x, sigma, theta, flag)})
  CovMat[2:(Nknots+1),1]<- t(CovMat[1,2:(Nknots+1)])
  CovMat[2:(Nknots+1),2:(Nknots+1)]<- mkpp(d[1:Nknots,1:Nknots], sigma, theta, flag)
  return(CovMat)
}


ESS_xi_tau = function(xi, xi_0, beta, tau, tau_0, logL, H_Z, 
                  h_y,X_new,test_knots){
  thetamin = 0
  thetamax = 2*pi
  u = runif(1)
  logy = logL + log(u)
  theta = runif(1, thetamin, thetamax)
  thetamin = theta - 2*pi
  thetamax = theta
  xistar = xi*cos(theta) + xi_0*sin(theta)
  taustar = tau*cos(theta) + tau_0*sin(theta)
  # xistar = xistar - t(kronecker(rep(1,knot_length),t(rowMeans(xistar)))) - 
  #   kronecker(rep(1,knot_length),t(colMeans(xistar))) + mean(xistar)
  logLstar = logcond_likelihood(xistar,beta,taustar,h_loc,X_new,h_y,test_knots,H_Z)$L
  
  while(logLstar<=logy){
    if(theta<0){ thetamin = theta }
    else{ thetamax = theta }
    theta = runif(1, thetamin, thetamax)
    xistar = xi*cos(theta) + xi_0*sin(theta)
    taustar = tau*cos(theta) + tau_0*sin(theta)
    # xistar = xistar - t(kronecker(rep(1,knot_length),t(rowMeans(xistar)))) - 
    #   kronecker(rep(1,knot_length),t(colMeans(xistar))) + mean(xistar)
    logLstar = logcond_likelihood(xistar,beta,taustar,h_loc,X_new,h_y,test_knots,H_Z)$L
  }
  return(list(xi = xistar, tau = taustar, logL = logLstar))       
}

# ESS_tau = function(xi, beta, tau,tau_0, logL, H_Z, 
#                   h_y,X_new,test_knots){
#   thetamin = 0
#   thetamax = 2*pi
#   u = runif(1)
#   logy = logL + log(u)
#   theta = runif(1, thetamin, thetamax)
#   thetamin = theta - 2*pi
#   thetamax = theta
#   taustar = tau*cos(theta) + tau_0*sin(theta)
#   #taustar = taustar - t(kronecker(rep(1,knot_length),t(rowMeans(taustar)))) - 
#   #  kronecker(rep(1,knot_length),t(colMeans(taustar))) + mean(taustar)
#   logLstar = logcond_likelihood(xi,beta,taustar,h_loc,X_new,h_y,test_knots,H_Z)$L
#   
#   while(logLstar<=logy){
#     if(theta<0){ thetamin = theta }
#     else{ thetamax = theta }
#     theta = runif(1, thetamin, thetamax)
#     taustar = tau*cos(theta) + tau_0*sin(theta)
#     # taustar = taustar - t(kronecker(rep(1,knot_length),t(rowMeans(taustar)))) - 
#     #   kronecker(rep(1,knot_length),t(colMeans(taustar))) + mean(taustar)
#     logLstar = logcond_likelihood(xi,beta,taustar,h_loc,X_new,h_y,test_knots,H_Z)$L
#   }
#   return(list(tau = taustar, logL = logLstar))       
# }

ESS_beta = function(xi, beta, beta_0,tau, logL, 
                    h_y,X_new,test_knots){
  thetamin = 0
  thetamax = 2*pi
  u = runif(1)
  logy = logL + log(u)
  theta = runif(1, thetamin, thetamax)
  thetamin = theta - 2*pi
  thetamax = theta
  betastar = beta*cos(theta) + beta_0*sin(theta)
  betastar = t(t(betastar)/sqrt(colSums(betastar * betastar)))
  bo<- logcond_likelihood_beta(xi,betastar,tau,h_loc,X_new,h_y,test_knots)
  logLstar = bo$L
  H_Z_mat = bo$H_Z.val
  
  while(logLstar<=logy){
    if(theta<0){ thetamin = theta }
    else{ thetamax = theta }
    theta = runif(1, thetamin, thetamax)
    betastar = beta*cos(theta) + beta_0*sin(theta)
    betastar = t(t(betastar)/sqrt(colSums(betastar * betastar)))
    bo<- logcond_likelihood_beta(xi,betastar,tau,h_loc,X_new,h_y,test_knots)
    logLstar = bo$L
  }
  H_Z_mat = bo$H_Z.val
  return(list(beta = betastar, logL = logLstar, H_Z = H_Z_mat))       
}

get_updates<- function(xi,beta,tau,h_loc,X,y,nmcmc = 2500,burnin = 2500,
                       thining = 5, test_knots,verbose=T){
  n <- nrow(X) #no. of individuals
  K <- ncol(y) #no. of Location
  p <- ncol(X) #no. of covariates
  
  knot_length = length(test_knots)
  XI_mat<- array(0, dim = c(knot_length,knot_length,nmcmc/thining))
  BETA_mat<- array(0,dim = c(p, K, nmcmc/thining))
  TAU_mat<- array(0, dim = c(knot_length,knot_length,nmcmc/thining))
  LOG_Like<- rep(0,nmcmc/thining)
  
  # max_norm<- max(sqrt(rowSums(X*X)))
  # X<- X/max_norm
  X_new<- kronecker(diag(K),X)
  
  y_new<- matrix(as.vector(y),nc=1)
  h_y<- hBasis(y_new,test_knots)
  intsize<- 1/(length(test_knots)-1)
  #intsize2<- 1/(length(test_knots2)-1)
  #test_knots1 = test_knots

  logL<- logcond_likelihood_beta(xi = xi,
                                 beta = beta, 
                                 tau = tau,
                                 h_loc=h_loc,
                                 X_new = X_new,
                                 h_y = h_y,
                                 test_knots = test_knots)$L
  
  t0 = proc.time()
  set.seed(728)
  for(i in 1:(nmcmc+burnin)){
    if(i%%100==0 & verbose==T) print(i)
    #####beta_update########
    beta_0<- matrix(as.vector(mvtnorm::rmvnorm(1,mean = rep(0, p*K))),nc=K)
    beta_0<- t(t(beta_0)/sqrt(colSums(beta_0 * beta_0)))
    result_beta <- ESS_beta(xi, beta, beta_0,tau, logL, 
                            h_y,X_new,test_knots)
    beta = result_beta$beta
    logL = result_beta$logL
    H_Z = result_beta$H_Z
    #####xi_tau_update########
    xi_0<- matrix(CHOL1%*%rnorm(knot_length^2),nc=knot_length,byrow = T)
    tau_0 =  matrix(CHOL1%*%rnorm(knot_length^2),nc=knot_length,byrow = T)
    # xi_0 = xi_0 - t(kronecker(rep(1,knot_length),t(rowMeans(xi_0)))) - 
    #   kronecker(rep(1,knot_length),t(colMeans(xi_0))) + mean(xi_0)
    result_xi_tau = ESS_xi_tau(xi, xi_0, beta,tau,tau_0, logL, H_Z, 
                       h_y,X_new,test_knots)
    xi = result_xi_tau$xi
    tau = result_xi_tau$tau
    logL = result_xi_tau$logL
    #####tau_update#######
    # tau_0 =  matrix(CHOL1%*%rnorm(knot_length^2),nc=knot_length,byrow = T)
    # # tau_0 = tau_0 - t(kronecker(rep(1,knot_length),t(rowMeans(tau_0)))) - 
    # #   kronecker(rep(1,knot_length),t(colMeans(tau_0))) + mean(tau_0)
    # result_tau = ESS_tau(xi, beta, tau, tau_0, logL, H_Z, 
    #                    h_y,X_new,test_knots)
    # tau = result_tau$tau
    # logL = result_tau$logL
    ######################
    #print(i)
    ######################
    if(i>burnin && i%%thining==0){
      XI_mat[  ,,(i-burnin)/thining]<- xi
      BETA_mat[  ,,(i-burnin)/thining]<- beta
      TAU_mat[  ,,(i-burnin)/thining] = tau
      LOG_Like[(i-burnin)/thining]<- logL
    }
  }
  print(proc.time()-t0)
  return(list(XI_mat.out = XI_mat, 
              BETA_mat.out = BETA_mat, 
              TAU_mat.out=TAU_mat,
              LOG_Like.out = LOG_Like))
}

####################################
##### Posterior Function ###########
####################################

PostFun_loc<- function(XI_mat, y_grid, X_new, BETA_mat, TAU_mat,
                       location,test_knots){
  intsize = 1/(length(test_knots)-1)
  if(ncol(BETA_mat) < location) 
  {return("Location is wrong!")
    break}
  K <- dim(BETA_mat)[2] #no. of Location
  p <- dim(BETA_mat)[1] #no. of covariates
  Niter<- dim(BETA_mat)[3]
  grid_point<- length(y_grid)
  samplesize = nrow(X_new)
  grid_mat<- array(0, dim = c(grid_point,Niter,samplesize))
  for(i in 1:Niter){
    xi<- XI_mat[,,i]
    beta<- as.matrix(BETA_mat[, location, i])
    tau<- TAU_mat[,,i]
    beta_tilde <- t(t(beta)/sqrt(colSums(beta * beta)))
    Z<- (X_new%*%beta_tilde + 1)/2
    xi_tilde = (xi%*%t(hBasis(Z,test_knots))) +
      kronecker(matrix(rep(1,samplesize),nr=1),
                tau%*%t(hBasis((location-1)/(K-1),test_knots)))
    top = hBasis(y_grid,test_knots)%*%xi_tilde
    xi_tilde.max = max(xi_tilde)
    A1 = xi_tilde[2:nrow(xi_tilde),] - xi_tilde.max
    A2 = xi_tilde[1:(nrow(xi_tilde)-1),] - xi_tilde.max
    const = log(colSums((exp(A2)-exp(A1))/((A2-A1)/intsize))) + xi_tilde.max
    grid_mat[,i,]<- exp(top - kronecker(rep(1,length(y_grid)),t(const)))
  }
  return(grid_mat)
}
######### Calculate Quantiles ################
get_quantile = function(xi, tau, beta, X_pred, test_knots, prob=0.5){
  knot_length = length(test_knots)
  intsize = 1/(knot_length-1)
  p=nrow(beta)
  K=ncol(beta)
  X_new = kronecker(diag(K),X_pred)
  n=nrow(X_new)/K
  location_grid<- seq(0,K-1,by=1)/(K-1)
  h_loc<-hBasis(location_grid,test_knots)
  beta_tilde <- t(t(beta)/sqrt(colSums(beta * beta)))
  beta_tilde_aug<- matrix(as.vector(beta_tilde),nc=1)
  Z<- (X_new%*%beta_tilde_aug + 1)/2
  delta = (xi%*%t(hBasis(Z,test_knots))) +
    kronecker(tau%*%t(h_loc),matrix(rep(1,n),nr=1))
  delta.max = max(delta)
  A2 = delta[2:nrow(delta),]
  A1 = delta[1:(nrow(delta)-1),]
  
  A3 = exp(delta.max)*(exp(A2-delta.max)-exp(A1-delta.max))/((A2-A1)/intsize)
  
  A4 = apply(A3,MARGIN = 2, function(x){cumsum(x)})
  
  constant = A4[nrow(A4),]
  
  #u = rep(prob,n*K)
  
  s=rep(0,n*K)
  for(i in 1:(n*K)){
    alpha = prob
    const = constant[i]
    j_0 = max(which(A4[,i]<alpha * const),0)+1
    if(j_0>1)
      s[i] = test_knots[j_0+1] + (log((alpha * const - A4[(j_0-1),i]) * 
                                        (((A2-A1)/intsize)[j_0,i]) + exp(delta[j_0,i])) - 
                                    delta[(j_0+1),i])/(((A2-A1)/intsize)[j_0,i])
    else
      s[i] = test_knots[j_0+1] + (log((alpha * const) * 
                                        (((A2-A1)/intsize)[j_0,i]) + exp(delta[j_0,i])) - 
                                    delta[(j_0+1),i])/(((A2-A1)/intsize)[j_0,i])
  }
  return(matrix(s,nc=K,byrow = F))
}
