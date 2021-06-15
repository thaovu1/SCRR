library(pracma)
library(MASS)
library(e1071)
library(caret)
library(glmnet)
library(doParallel)
#set directory to load in source functions
source("~/source_functions.R")


set.seed(0807)
s = 0.0005
s2 = 0.002
k = 50
e  = 0.0001
s1 = 0.001
n_extra = 190


ppm = seq(0.9, 9.2, by = 0.001)
ppm = round(ppm, 3)
n = length(ppm)
u1 = which(ppm == 1.874)
u2 = which(ppm == 2.1)
# peak locations and intensity
p1 = 1.924
p3 =  2
p4 = 2.05
p5 = 1.8
p6 = c(1.708, 1.828)
p7 = c(1.727, 1.728, 1.822)
p8 = c(1.756, 1.789, 1.790, 1.791)
p9 = c(1.827, 1.956, 2.077)
p10 = c(1.801, 1.802, 1.803, 1.804, 1.878, 1.879)
p11 = c(1.977, 1.978)
# intensity
i6 = c(1, 0.6)
i7 = c(0.428, 0.762, 1)
i8 = c(0.6, 0.7,1, 0.85)
i9 = c(1, 0.5, 0.6)
i10 = c(0.55, 0.9, 0.8, 0.4, 1,0.7)
i11 = c(0.7, 1)
test_peak = list(p1,p3,p4,p5, p6,p7,p8,p9,p10,p11)
test_peak = lapply(test_peak, function(x) round(x, 3))
int_z = list(1,1/3,1/2,1,i6,i7,i8,i9,i10,i11)
test_peak_id = lapply(test_peak, function(x)(sapply(x, function(x) which(ppm == x))))
peaks_id1 = sort(unlist(test_peak_id))
peaks_id2 = id_create(peaks_id1, 10 , n)
c1 = which(ppm == 4.6)
c2 = which(ppm == 4.821)
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
np_total = sum(np_all)
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
    refset[,z] = sim_cauchy(loc_list1[[z]], int_list1[[z]], ppm, s2)
  }
}
shape_parm = c(rep(list(s),length(test_peak)), rep(list(s2),n_extra))
# peak indexes
peaks_id = lapply(loc_list1, function(x)(sapply(x, function(x) which(ppm == x))))
peaks_idx1 = sort(unlist(peaks_id))

# REPLICATE: Proposed method, LASSO, Adaptive Lasso, Elastic net
set.seed(0807)
num_iters = 5
n_repeat = 10 #200
sig = 1/3 #sigma_0
J1 = ncol(refset)
s_r = 10 #shift_range
s_amt = c() #shift_amount
t_v1 = t_v2 = matrix(0, n_repeat, J1)
start.time = Sys.time()
for (v in 1:n_repeat){
  if(v %% 50==0){
    cat(paste0("iteration: ", v, "\n"))
  }
  d1 = sample(c(-s_r:s_r), 1)
  d = abs(d1)+1
  shift = d1*0.001
  # shifted spectrum 1
  s1_shifted = sim_cauchy(round(loc_list1[[1]]+shift,3), int_list1[[1]], ppm, s)
  true_sim = s1_shifted + refset[,2] + refset[,3] 
  noise_std = trapz(ppm, true_sim)
  noise = rnorm(n, 0, noise_std)
  real_sim = true_sim + noise
  real_sim = as.matrix(real_sim)
  #threshold
  peak_use = which(real_sim > 1.07*max(real_sim[c2:which.max(ppm)]))
  # search through reference library to remove metabolites that don't have at least one peak in peak_use
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
 
  beta_opt = res_opt$beta_hat
  ref_idx = which(beta_opt > 0)
  # SHIFTING CORRECTION
  if (length(ref_idx) < 1){
    next
  }
  w_opt = res_opt$weight_matrix
  w_max = apply(w_opt, 1, max)
  w_id = which(w_max > 0.05)
  w_rm = which(w_opt[w_id,d+1] > 0.2)
  if (length(w_rm) > 0) {w_id1 <- w_id[-w_rm]
  } else {w_id1 <- w_id}
  
  # only consider i's that don't have max weight at the center
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
J1 = ncol(refset)
s_r = 10 #shift_range

t_v3 = t_v4 = matrix(0, n_repeat, J1)
start.time = Sys.time()
for (v in 1:n_repeat){
  if(v %% 50==0){
    cat(paste0("iteration: ", v, "\n"))
  }
  d1 = sample(c(-s_r:s_r), 1)
  d = abs(d1)+1
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

