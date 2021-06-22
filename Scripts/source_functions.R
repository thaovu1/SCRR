# SOURCE FUNCTIONS
library(MASS)
library(pracma)


# Simulating NMR spectra
sim_cauchy = function(x, y, ppm, s){
  #x: location of peaks
  #y: corresponding intensity
  #ppm: input vector of chemical shift
  #s: shape parameter
  p = length(ppm)
  l = length(x)
  index = c()
  m = matrix(0, l, p)
  if(length(y) >1){
    y = y/max(y)}
  for (i in 1:l){
    m[i,] = dcauchy(ppm, x[i], s)/(s*pi)
    index[i] = which(ppm == x[i])
  }
  coef = t(m[,index])
  solution = as.numeric(t(ginv(coef) %*% y))
  s1 = crossprod(m,solution)
  return (s1) 
}

# create important region
id_create = function(id, d = 10, n){
  #n: total number of data points
  res = c()
  np = length(id)
  for (k in 1:np){
    i = id[k]
    pwindow = seq(max(1, i-d), min(n, i+d))
    res = sort(union(res, pwindow))
  }
  return(res)
}

# exponential part used in weight function
p_s1 = function(x,y, sig){
  res = exp((-x^2 + y^2)/(2*sig^2))
  return(res)
}

# weight function used in the inner summation
w_weight = function(s, sig = 0.01){
  # x1 should be a vector
  k = length(s)
  w = c()
  for (i in 1:k){
    w[i] = 1/sum(p_s1(s, s[i], sig))
  }
  return(as.matrix(w))
}

# soft-threshold operator with non-negativity constraint
soft_threshold = function(rho, lambda, z){
  res = 0
  if (rho < -lambda){
    res = max(0, (rho + lambda)/z)
  }  else if (rho > lambda){
    res = max(0,(rho - lambda)/z)
  }
  return (res)
}

# "1/N" variance form used to standardize X
mysd <- function(y) sqrt(sum((y-mean(y))^2)/length(y))

# lambda_max after thresholding
l_thres = function(x, y, peaks_id, d = 10){
  #x: reference library matrix
  #y: response spectrum
  #peaks_id: indexes of signals after thresholding
  #d: search window size
  n = nrow(x)
  x = scale(x,scale = apply(x,2,mysd))
  p = ncol(x)
  y = scale(y, scale = F)
  h = peaks_id
  n1 = length(h)
  result = c()
  for (j in 1:p){
    total = 0
    for (i in h){
      idx = seq(max(1, i-d), min(n, i+d))
      nrep = length(idx)
      x_j = x[idx, j]
      y_obs = rep(y[i], nrep)
      total = total + (1/nrep)*crossprod(y_obs, x_j)
    }
    result[j] = total
  }
  max_l = max(abs(result))/n1
  return(max_l)
}

#find neighbor indexes
idx_seq = function(id, d, n){
  m = length(id)
  # length 
  l = 2*d +1
  res = matrix(0, l, m)
  for (i in 1:m){
    s = seq(max(1, id[i]-d), min(n, id[i]+d))
    res[1:length(s),i] = s
  }
  return(res)
}

#range of peak with center at i
idx_seq1 = function(id, w, n){
  res = sapply(id, function(x) seq(x-w,x+w))
  return(res)
}


#regression parameter estimation
thres_coordinate = function(theta, X, y, peaks_id, lambda = 0.01, num_iters = 100, d = 10, 
                  sig_i, sc = T){ 
  #theta: initial parameter values
  #X: matrix of reference spectra
  #y: vector of response spectrum
  #peaks_id: indexes for signals after thresholding
  #lambda: penalty value
  #num_iters: number of iterations to run coordinate descent
  #d: search window size
  #sig_i: value of weight distributor sigma_0
  #sc: indicator for scaling
  n = nrow(X)
  if (sc == T){
    X = scale(X, scale = apply(X, 2, mysd))
    y = scale(y, scale = F)
    }  else { 
    X = X 
    y = y} 
  J = ncol(X)
  h = peaks_id
  n1 = length(peaks_id)
  for (v in 1:num_iters){
    for (j in 1:J){
      rho = 0
      z = 0
      res_w = matrix(0, n, 2*d+1)
      for (i in h){
        idx = seq(max(1, i-d), min(n, i+d))
        nrep = length(idx)
        X_j = X[idx, j]
        y_pred = X[idx,] %*% theta
        y_obs = rep(y[i], nrep)
        resi = y_obs - y_pred
        rss = resi^2
        w = w_weight(rss, sig_i)
        res_w[i,1:nrep] = w
        rho = rho + crossprod(w*X_j, y_obs - y_pred + theta[j]*X_j)
        z = z + crossprod(w, X_j^2)
      }
      theta[j] = soft_threshold(rho/n1, lambda, z/n1)
      }
  }
  return(list("beta_hat" = theta, "weight_matrix" = res_w))
}

# Cross validation
cv_thres_coordinate = function(theta, X, y, peaks_id, lambda = 0.01, num_iters = 100, d = 10, 
                    sig_i, nfold = 10){ 
  n = nrow(X)
  X = scale(X, scale = apply(X, 2, mysd))
  y = scale(y, scale = F)
  J = ncol(X)
  np = length(peaks_id)
  fold <- cut(seq(1,np),breaks = nfold,labels=FALSE)
  newfold = sample(fold)
  # prediction error
  pe = c()
  for (u in 1:nfold){
    testid = peaks_id[which(newfold == u,arr.ind=TRUE)]
    trainid = peaks_id[-which(newfold == u,arr.ind=TRUE)]
    test_X = X[testid,]
    test_y = y[testid]
    ntest = nrow(test_X)
    temp = thres_coordinate(theta, X, y, trainid,lambda, num_iters,
                            d, sig_i, sc = F)
    theta_hat = temp$beta_hat
    rss = crossprod(test_y - test_X %*% theta_hat, 
                    test_y - test_X %*% theta_hat)
    w_hat = w_weight(rss, sig_i)
    # loss function according to Eq. (2.6)
    loss = crossprod(w_hat, rss)
    pe[u] = loss
  }
  mpe = mean(pe)
  return(mpe)
}



 
