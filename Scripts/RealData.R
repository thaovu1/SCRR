rm(list = ls())
library(pracma)
library(MASS)
library(e1071)
library(caret)
library(glmnet)
library(doParallel)
# specify the path which has source code file to load functions
source("source_functions.R")
# specify the path to Input directory to load data and reference libraries
# load chemical shifts
load("ppm.RData")
# load Experimental mixture 
load("Exp_data.RData")
# true metabolite indexes
load("Exp_truth.RData")
# load Hart et al. serum samples
load("Hart_data.RData")
# true metabolite indexes
load("Hart_truth.RData")
# load all the libraries
load("spec_lib61.RData")
load("spec_lib101.RData")
load("spec_lib104.RData")
load("spec_lib200.RData")
# extract reference matrices and corresponding peak locations and intensities
refset61 = spec_lib61$refset61
peaks_loc61 = spec_lib61$peaks_loc61
peaks_int101 = spec_lib61$peaks_int61

refset101 = spec_lib101$refset101
peaks_loc101 = spec_lib101$peaks_loc101
peaks_int101 = spec_lib101$peaks_int101

refset104 = spec_lib104$refset104
peaks_loc104 = spec_lib104$peaks_loc104
peaks_int104 = spec_lib104$peaks_int104

refset200 = spec_lib200$refset200
peaks_loc200 = spec_lib200$peaks_loc200
peaks_int200 = spec_lib200$peaks_int200


#=========================================================================================#
# with library of 61
#=========================================================================================#
# change these accordingly based on which library being considered
refset = refset61
peaks_loc = peaks_loc61
peaks_int = peaks_int61
# mixture 1
mix1 = as.numeric(Exp_data[1,])
mix_truth = Exp_truth$Mix1

#predefined parameters
d = 10  #search window size
k = 50  #number of lambda values to evaluate
e  = 0.0001 #arbitrary number used to create a sequence of lambda values
sig = 1/3  #sigma_0: weight distributor parameter
s = 0.002 # peak shape parameter
num_iters = 5  #number of iterations
thres = 0.07*trapz(ppm, mix1) 
peak_use = which(mix1 > thres)

l_max = l_thres(refset, mix1, peak_use, d)
l_seq = exp(seq(log(e*l_max), log(l_max), length.out = k))
l_seq = rev(l_seq)

J  = ncol(refset) 
beta_init =  matrix(0,J,1)
#=========================================================================================#
# Cross-validation
#=========================================================================================#
# set.seed()
mpe = c()
for (a in 1:k){
  mpe[a] = cv_thres_coordinate(beta_init, refset, mix1, peak_use, lambda = l_seq[a],
                               num_iters, d, sig, nfold = 5)
}
lambda_opt = l_seq[which.min(mpe)]

res_opt = thres_coordinate(beta_init, refset, mix1, peak_use, lambda = lambda_opt,
                           num_iters, d, sig, sc = T)

beta_opt = res_opt$beta_hat
ref_idx = which(beta_opt > 0)

# SHIFTING CORRECTION
w_opt = res_opt$weight_matrix
w_max = apply(w_opt, 1, max)
w_id = which(w_max > 0.05)
# only look at i's that don't have max weight at the center
w_id1 = w_id
if (length(w_id) >1){
  w_max_id = apply(w_opt[w_id1,] ,1, which.max)
}else{w_max_id = which.max(w_opt[w_id1,])}

w_max_id1 = w_id1[which(w_max_id != d+1)]
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
new_loc = peaks_loc[ref_idx]
new_int = peaks_int[ref_idx]

new_ref = as.matrix(refset[,ref_idx])
new_loc2 = new_loc
for (i in 1:length(new_loc2)){
  result = lapply(new_loc2[[i]], function(x) which(h_id == x))
  if (sum(lengths(result))>0){
    m = which(lengths(result) >0)
    for (g in m){
      new_loc2[[i]][g] = peak_change[result[[g]]][1]
    }
    
    new_ref[,i] = sim_cauchy(new_loc2[[i]], unlist(new_int[[i]]), ppm, s)
  }
}
# SECOND-STAGE FITTING
mdl = lm(formula = mix1 ~ new_ref - 1 )
beta_tilde = mdl$coefficients

#calculate TP, FP, FN

FP = length(setdiff(ref_idx[which(beta_tilde>0)], mix_truth))
FN = length(setdiff(mix_truth, ref_idx[which(beta_tilde>0)]))
TP = length(intersect(mix_truth, ref_idx[which(beta_tilde>0)]))
TN = J - (FP + FN + TP)



#=========================================================================================#
# Hart et al. samples
#=========================================================================================#
refset = refset104
peaks_loc = peaks_loc104
peaks_int = peaks_int104
# sample 1
mix1 = as.numeric(Hart_data[1,])
mix_truth = Hart_truth
# continue like above
#predefined parameters
d = 10  #search window size
k = 50  #number of lambda values to evaluate
e  = 0.0001 #arbitrary number used to create a sequence of lambda values
sig = 1/3  #sigma_0: weight distributor parameter
s = 0.002 # peak shape parameter
num_iters = 5  #number of iterations
thres = 0.07*trapz(ppm, mix1) 
peak_use = which(mix1 > thres)

l_max = l_thres(refset, mix1, peak_use, d)
l_seq = exp(seq(log(e*l_max), log(l_max), length.out = k))
l_seq = rev(l_seq)

J  = ncol(refset) 
beta_init =  matrix(0,J,1)
#=========================================================================================#
# Cross-validation
#=========================================================================================#
# set.seed()
mpe = c()
for (a in 1:k){
  mpe[a] = cv_thres_coordinate(beta_init, refset, mix1, peak_use, lambda = l_seq[a],
                               num_iters, d, sig, nfold = 5)
}
lambda_opt = l_seq[which.min(mpe)]

res_opt = thres_coordinate(beta_init, refset, mix1, peak_use, lambda = lambda_opt,
                           num_iters, d, sig, sc = T)

beta_opt = res_opt$beta_hat
ref_idx = which(beta_opt > 0)

# SHIFTING CORRECTION
w_opt = res_opt$weight_matrix
w_max = apply(w_opt, 1, max)
w_id = which(w_max > 0.05)
# only look at i's that don't have max weight at the center
w_id1 = w_id
if (length(w_id) >1){
  w_max_id = apply(w_opt[w_id1,] ,1, which.max)
}else{w_max_id = which.max(w_opt[w_id1,])}

w_max_id1 = w_id1[which(w_max_id != d+1)]
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
new_loc = peaks_loc[ref_idx]
new_int = peaks_int[ref_idx]

new_ref = as.matrix(refset[,ref_idx])
new_loc2 = new_loc
for (i in 1:length(new_loc2)){
  result = lapply(new_loc2[[i]], function(x) which(h_id == x))
  if (sum(lengths(result))>0){
    m = which(lengths(result) >0)
    for (g in m){
      new_loc2[[i]][g] = peak_change[result[[g]]][1]
    }
    
    new_ref[,i] = sim_cauchy(new_loc2[[i]], unlist(new_int[[i]]), ppm, s)
  }
}
# SECOND-STAGE FITTING
mdl = lm(formula = mix1 ~ new_ref - 1 )
beta_tilde = mdl$coefficients

#calculate TP, FP, FN

FP = length(setdiff(ref_idx[which(beta_tilde>0)], mix_truth))
FN = length(setdiff(mix_truth, ref_idx[which(beta_tilde>0)]))
TP = length(intersect(mix_truth, ref_idx[which(beta_tilde>0)]))
TN = J - (FP + FN + TP)




