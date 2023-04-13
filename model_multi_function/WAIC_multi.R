den.eval_model.2<- function(XI_mat, y, X, BETA_mat,
                            location,test_knots){
  knot_length = length(test_knots)
  intsize = 1/(knot_length-1)
  if(ncol(BETA_mat) < location) 
  {return("Location is wrong!")
    break}
  K <- dim(BETA_mat)[2] #no. of Location
  p <- dim(BETA_mat)[1] #no. of covariates
  Niter<- dim(BETA_mat)[3]
  grid_point<- length(y)
  samplesize = nrow(X)
  cond.den<- array(0, dim = c(samplesize,Niter))
  for(i in 1:Niter){
    xi<- XI_mat[(knot_length*(location-1)+1):(knot_length*location),,i]
    beta<- as.matrix(BETA_mat[, location, i])
    beta_tilde <- t(t(beta)/sqrt(colSums(beta * beta)))
    Z<- (X%*%beta_tilde + 1)/2
    xi_tilde = (xi%*%t(hBasis(Z,test_knots)))
    top = hBasis(y,test_knots)%*%xi_tilde
    xi_tilde.max = max(xi_tilde)
    A1 = xi_tilde[2:nrow(xi_tilde),] - xi_tilde.max
    A2 = xi_tilde[1:(nrow(xi_tilde)-1),] - xi_tilde.max
    const = log(colSums((exp(A2)-exp(A1))/((A2-A1)/intsize))) + xi_tilde.max
    cond.den[,i]<- diag(exp(top - kronecker(rep(1,length(y)),t(const))))
  }
  return(cond.den)
}

waic.eval_model.2 = function(XI_mat, y, X, BETA_mat, test_knots){
  K = ncol(y)
  samplesize = nrow(X)
  Niter = dim(BETA_mat)[3]
  dummy.mat = data.frame()
  for(i in 1:K){
    A = as.data.frame(den.eval_model.2(XI_mat, y[,i], X, BETA_mat,
                                       location=i,test_knots))
    dummy.mat = rbind(dummy.mat,A)
  }
  dummy.mat = as.matrix(dummy.mat)
  pwaic1 = 2*sum(log(rowMeans(dummy.mat)) - rowMeans(log(dummy.mat)))
  pwaic2 = sum(apply(log(dummy.mat),1,var))
  lppd = sum(log(apply(dummy.mat,1,mean)))
  return(list( elppd1 = -2*(lppd-pwaic1), elppd2 = -2*(lppd-pwaic2)))
}
