library(pracma)
library(MASS)
library(e1071)
library(caret)
library(glmnet)
library(doParallel)
source("~/Documents/MVAPACK related/Shifting error regularization model paper /SRCC_code/source_functions.R")

set.seed(0807)
s2 = 0.002
#s = 0.001
s = 0.0005
#d = 10
k = 50
e  = 0.0001
s1 = 0.001
n_extra = 190
# ppm
#(x<0.8|x>4.599 & x<4.821|x>9.499)
ppm = seq(0.9, 9.2, by = 0.001)
ppm = round(ppm, 3)
n = length(ppm)

peak1 = c(1.8, 1.803, 1.87)
peak2 = c(1.759, 1.764, 1.951,1.966,1.980,1.995)
peak3 = c(1.903, 2.03)
peak4 = c(1.821, 1.841, 1.947, 1.977)
peak5 = c(1.925, 1.935,1.94, 2.00, 2.011,2.021)
peak6 = c(1.705, 1.715)
peak7= c(1.645, 1.646, 1.647, 2.877)
peak8 = c(2.271, 2.272, 2.273, 2.274)
peak9=c(1.583,2.647)
peak10 = 2.513
# intensities
int1 = c(0.8, 0.9, 1)
#int2 = c(0.9806, 1.0000, 0.0631, 0.2025, 0.2030, 0.0649)
int2 = c(0.7, 1.0000, 0.0631, 0.2025, 0.2030, 0.0649)
int3 = c(1, 0.4374)
int4 = c(0.5785,1,0.9238,0.6328)
int5 = c(0.9128,1,0.5884,0.9641,0.7559,0.8481)
int6 = c(0.3, 1)
int7= c(0.557, 0.965, 0.446, 1)
int8 = c(0.235, 0.672, 0.453, 1)
int9 = c(0.37, 1)
int10 = 1
# combine peaks
test_peak = list(peak1,peak2, peak3, peak4, peak5, peak6,peak7,peak8,peak9,peak10)
test_peak = lapply(test_peak, function(x) round(x, 3))
int_z = list(int1,int2,int3,int4,int5,int6,int7,int8,int9,int10)
test_peak_id  = lapply(test_peak, function(x)(sapply(x, function(x) which(ppm == x))))
peaks_id1 = sort(unlist(test_peak_id))
# create a region such that all other ref spectra within +_ ppm of above peaks
peaks_id2 = id_create(peaks_id1, 10 , n)
# dont use first and last 10 ppm to generate peaks
c1 = which(ppm == 4.6)
c2 = which(ppm == 4.821)
#ppm1 = ppm[-c(head(1:n,100),tail(1:n,100), u1:u2, c1:c2)] 
ppm1 = ppm[-c(head(1:n,100),tail(1:n,100), c1:c2, peaks_id2)] 


loc_list = list() 
int_list = list()
np = list()

for (i in 1:n_extra){
  # number of peaks per metabolite
  np[[i]] = sample(c(1:3),1)
}
# total # peaks
np_all = unlist(np)
index_lower = cumsum(np_all) - np_all + 1
index_upper = cumsum(np_all)
index_mat = cbind(index_lower, index_upper)
np_total = sum(unlist(np))
# all peak locations
all_peaks = sample(ppm1, np_total, replace = F) #to prevent duplication
for (i in 1:n_extra){
  l = np[[i]]
  loc_list[[i]] = all_peaks[index_mat[i,1]:index_mat[i,2]]
  int_list[[i]] = round(runif(l, 0.1, 1),4)
}

loc_list1 = c(test_peak, loc_list)
int_list1 = c(int_z, int_list)

ntotal = length(loc_list1)
refset = matrix(0, n, ntotal)
for (z in 1:ntotal){
  if(z < (length(test_peak)+1)){
    refset[,z] = sim_cauchy(loc_list1[[z]], int_list1[[z]], ppm, s)
  }
  else{
    # s2 or s?
    refset[,z] = sim_cauchy(loc_list1[[z]], int_list1[[z]], ppm, s)
  }
}
shape_parm = c(rep(list(s),length(test_peak)), rep(list(s),n_extra))
# peak indexes
peaks_id = lapply(loc_list1, function(x)(sapply(x, function(x) which(ppm == x))))
peaks_idx1 = sort(unlist(peaks_id))
# repeat
set.seed(0807)
n_repeat = 200
sig = 1/3
s_r = 40 #shift range
J1 = ncol(refset)
t_v1 = t_v2 = matrix(0, n_repeat, J1)
start.time = Sys.time()
for (v in 1:n_repeat){
  if(v %% 50==0){
    cat(paste0("iteration: ", v, "\n"))
  }
  d1 = sample(c(-s_r:s_r), 1)
  d = s_r
  shift = d1*0.001
  peak2_shifted = loc_list1[[2]]
  peak2_shifted[1:2]= round(peak2_shifted[1:2] + shift, 3)
  s2_shifted = sim_cauchy(peak2_shifted, int_list1[[2]], ppm, s)
  true_sim = refset[,1] + s2_shifted +  refset[,3]  + refset[,4] + refset[,5]
  noise = rnorm(n, 0, 0.01)
  real_sim = true_sim + noise
  real_sim = as.matrix(real_sim)
  peak_use = which(real_sim > 1.07*max(real_sim[c2:which.max(ppm)]))
  peak_check = id_create(peak_use, d, n)
  check_id = list()
  for (i in 1:ntotal){
    check_id[[i]] = lapply(peaks_id[[i]], function(x) which(peak_check == x))
  }  
  comp_mix = lapply(check_id, function(x) sum(lengths(x)))
  comp_mix1 = which(unlist(comp_mix) > 0)
  ref_filter = refset[,comp_mix1]
  J = ncol(ref_filter)
  #initial values
  beta_init = matrix(0, J, 1)
  
  l_max = l_thres(ref_filter, real_sim, peak_use, d)
  l_seq = exp(seq(log(e*l_max), log(l_max), length.out = k))
  l_seq = rev(l_seq)
  
  # cross validation
  mpe = c()
  for (a in 1:k){
    mpe[a] = cv_thres_coordinate(beta_init, ref_filter, real_sim, peak_use, lambda = l_seq[a],
                                 num_iters, d, sig, nfold = 5)
  }
  lambda_opt = l_seq[which.min(mpe)]
  
  res_opt = thres_coordinate(beta_init, ref_filter, real_sim, peak_use, lambda = lambda_opt,
                             num_iters, d, sig, sc = T)

  # SHIFTING CORRECTION
  beta_opt = res_opt$beta_hat
  ref_idx = which(beta_opt > 0)
  if (length(ref_idx) < 1){
    next
  }
  w_opt = res_opt$weight_matrix
  w_max = apply(w_opt, 1, max)
  w_id = which(w_max > 0.05)
  w_rm = which(w_opt[w_id,d+1] > 0.2)
  if (length(w_rm) > 0) {w_id1 <- w_id[-w_rm]
  } else {w_id1 <- w_id}
  
  # only look at i's that don't have max weight at the center
  if (length(w_id) >1){
    w_max_id = apply(w_opt[w_id1,] ,1, which.max)
  }else{w_max_id = which.max(w_opt[w_id1,])}
  
  w_max_id1 = w_id1[which(w_max_id != d+1)]
  if (length(w_max_id1) < 1){
    mdl_ols = lm(formula = real_sim ~ ref_filter[,ref_idx] - 1 )
    beta_ols = mdl_ols$coefficients
    ref_idx1 = comp_mix1[ref_idx]
    t_v1[v, ref_idx1] = round(beta_ols,3)
    next
  }
  peak_change = ppm[w_max_id1]
  # neighbor indexes, (2d+1) x npeak
  idx = idx_seq(w_max_id1, d, n)
  idx1 = idx
  for (b in 1:ncol(idx)){
    # find nonzero elements to be replaced
    r = which(idx[,b] > 0)
    idx1[r,b] = ppm[idx[,b]]
  }
  # identify the shift
  if (length(w_max_id1) > 1){
    max_id = apply(w_opt[w_max_id1,] ,1, which.max)
  }else {max_id = which.max(w_opt[w_max_id1,])}
  # this should be index of true peak
  h_id = c()
  for (i in 1:length(w_max_id1)){
    h_id[i] = idx1[max_id[i],i]
  }
  new_loc = loc_list1[ref_idx]
  new_int = int_list1[ref_idx]
  #---------------------------------------------------------------------------
  new_ref = as.matrix(ref_filter[,ref_idx])
  new_loc2 = new_loc
  for (i in 1:length(new_loc2)){
    result = lapply(new_loc2[[i]], function(x) which(h_id == x))
    if (sum(lengths(result))>0){
      m = which(lengths(result) >0)
      for (g in m){
        new_loc2[[i]][g] = peak_change[result[[g]]][1]
      }
      new_ref[,i] = sim_cauchy(new_loc2[[i]], new_int[[i]], ppm, shape_parm[[ref_idx[i]]])
    }
  }
  
  # SECOND-STAGE FITTING
  
  mdl_ols = lm(formula = real_sim ~ new_ref - 1 )
  beta_ols = mdl_ols$coefficients
  ref_idx1 = comp_mix1[ref_idx]#comp_mix1[which(beta_ols)>0]
  t_v1[v, ref_idx1] = round(beta_ols,3)
  
  #--------------------------------------------------#
  # LASSO
  #--------------------------------------------------#
  cv2 = cv.glmnet(refset, real_sim, lambda = NULL, 
                  type.measure = "mse", alpha = 1, intercept = T, nlambda = 50, nfolds = 5)
  mdl2 = glmnet(refset, real_sim, lambda = cv2$lambda.min, intercept = T, lower.limits = 0)
  beta2 = predict(mdl2, type = "coefficients", s = cv2$lambda.min)
  beta2 = beta2[-1]
  #length(which(beta2>0))
  t_v2[v,] = beta2
}


end.time = Sys.time()
end.time - start.time

#--------------------------------------------------##--------------------------------------------------#
# ELASTIC NET AND ADAPTIVE LASSO
#--------------------------------------------------##--------------------------------------------------#
set.seed(0807)
n_repeat = 200
sig = 0.3
J1 = ncol(refset)
s_r = 10 #shift_range: 10, 20, 30, 40

t_v3 = t_v4 = matrix(0, n_repeat, J1)
start.time = Sys.time()
for (v in 1:n_repeat){
  if(v %% 50==0){
    cat(paste0("iteration: ", v, "\n"))
  }
  d1 = sample(c(-s_r:s_r), 1)
  d = s_r
  shift = d1*0.001
  # shifted spectrum 1
  s1_shifted = sim_cauchy(round(loc_list1[[1]]+shift,3), int_list1[[1]], ppm, s)
  true_sim = s1_shifted + refset[,2] + refset[,3] 
  noise_std = trapz(ppm, true_sim)
  noise = rnorm(n, 0, noise_std)
  real_sim = true_sim + noise
  real_sim = as.matrix(real_sim)
  #--------------------------------------------------#
  # ELASTIC NET 
  #--------------------------------------------------#
  a1 <- seq(0.1, 0.9, 0.025)
  search <- foreach(i = a1, .combine = rbind) %dopar% {
    cv <- cv.glmnet(refset, real_sim, type.measure = "mse", paralle = TRUE, alpha = i, 
                    intercept = T, nlambda = 50, nfolds = 5, lower.limits = 0)
    data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.min], lambda.min = cv$lambda.min, alpha = i)
  }
  cv3 <- search[search$cvm == min(search$cvm), ]
  mdl3 <- glmnet(refset, real_sim, lambda = cv3$lambda.min, alpha = cv3$alpha,intercept = T, lower.limits = 0)
  beta3 = predict(mdl3, type = "coefficients", s = cv3$lambda.min)
  beta3 = beta3[-1]
  t_v3[v,] = beta3
  #--------------------------------------------------#
  #ADAPTIVE LASSO
  #--------------------------------------------------#
  cv2 = cv.glmnet(refset, real_sim, lambda = NULL, 
                  type.measure = "mse", alpha = 1, intercept = T, nlambda = 50, nfolds = 5)
  mdl2 = glmnet(refset, real_sim, lambda = cv2$lambda.min, intercept = T, lower.limits = 0)
  beta2 = predict(mdl2, type = "coefficients", s = cv2$lambda.min)
  beta2 = beta2[-1]
  tau=1
  penalty.factor=abs(beta2+1/sqrt(n)^(-tau))
  cv_adapt = cv.glmnet(refset, real_sim,lambda = NULL, type.measure = "mse", alpha = 1,
                       intercept = T, nlambda = 50, nfolds = 5,penalty.factor = penalty.factor)
  mdl_adapt = glmnet(refset, real_sim, lambda = cv_adapt$lambda.min, intercept = T, lower.limits = 0)
  beta_adapt = predict(mdl_adapt, type = "coefficients", s = cv_adapt$lambda.min)
  beta_adapt = beta_adapt[-1]
  t_v4[v,] = beta_adapt
}
end.time = Sys.time()
end.time - start.time

#--------------------------------------------------------------------------------------------------------#
# ACCURACY, SENSITIVITY, SPECIFICITY
#--------------------------------------------------------------------------------------------------------#
# FP and FN
fp_v1 = apply(t_v1[,4:200], 1, function(x) length(which(x>0)))
fp_v2 = apply(t_v2[,4:200], 1, function(x) length(which(x>0)))
fp_v3 = apply(t_v3[,4:200], 1, function(x) length(which(x>0)))
fp_v4 = apply(t_v4[,4:200], 1, function(x) length(which(x>0)))

fn_v1 = apply(t_v1[,1:3], 1, function(x) length(which(x == 0)))
fn_v2 = apply(t_v2[,1:3], 1, function(x) length(which(x == 0)))
fn_v3 = apply(t_v3[,1:3], 1, function(x) length(which(x == 0)))
fn_v4 = apply(t_v4[,1:3], 1, function(x) length(which(x == 0)))

# TP and TN
tp_v1 = apply(t_v1[,1:3], 1, function(x) length(which(x>0)))
tp_v2 = apply(t_v2[,1:3], 1, function(x) length(which(x>0)))
tp_v3 = apply(t_v3[,1:3], 1, function(x) length(which(x>0)))
tp_v4 = apply(t_v4[,1:3], 1, function(x) length(which(x>0)))

tn_v1 = apply(t_v1[,4:200], 1, function(x) length(which(x==0)))
tn_v2 = apply(t_v2[,4:200], 1, function(x) length(which(x==0)))
tn_v3 = apply(t_v3[,4:200], 1, function(x) length(which(x==0)))
tn_v4 = apply(t_v4[,4:200], 1, function(x) length(which(x==0)))
#accuracy
A1 = (tp_v1 + tn_v1)/200
A2 = (tp_v2 + tn_v2)/200
A3 = (tp_v3 + tn_v3)/200
A4 = (tp_v4 + tn_v4)/200
#sensitivity
sen1 = tp_v1/(tp_v1+fn_v1)
sen2 =  tp_v2/(tp_v2+fn_v2)
sen3 =  tp_v3/(tp_v3+fn_v3)
sen4 =  tp_v4/(tp_v4+fn_v4)
# specificity
spec1 = tn_v1/(fp_v1+tn_v1)
spec2 = tn_v2/(fp_v2+tn_v2)
spec3 = tn_v3/(fp_v3+tn_v3)
spec4 = tn_v4/(fp_v4+tn_v4)

metric = data.frame("Accuracy" = round(c(mean(A1), mean(A2), mean(A3), mean(A4)),3),
                    "Sen" = round(c(mean(sen1), mean(sen2), mean(sen3), mean(sen4)), 3),
                    "Spec" = round(c(mean(spec1), mean(spec2), mean(spec3), mean(spec4)),3),
                    "sd_A" = round(c(sd(A1), sd(A2), sd(A3), sd(A4)),3),
                    "sd_sen" = round(c(sd(sen1), sd(sen2), sd(sen3),sd(sen4)),3),
                    "sd_spec" = round(c(sd(spec1), sd(spec2), sd(spec3), sd(spec4)),3))

#beta_tilde
final_beta = colMeans(t_v1[,1:3])
final_beta_sd = apply(t_v1[,1:3], 2, sd)

