library(MFPCA);library(tidyverse);library(MASS);require(zipfR);require(fOptions)
require(Matrix);library(parallel);library(foreach);library(doMC);library(rms);library(pec)
library(fda);library(survAUC);library(fda.usc)

source("ancillary.R")

#------- one run of the simulation under simulation with moderate censoring ------#

method = 'FPLS' # options: FPLS, sFPCA

n = 400;nnn = 300
ncr = 32
npc = 1 # number of basis to pick
tune_picked = 0.01 # tunning parameter for sFPCA, pre-specified in this example
lam_picked = 0.01 # tunning parameter for FPLS, pre-specified in this example
train_dat.id = tdat.id = read_rds('train_dat')
test_dat.id = read_rds('test_dat')

#----- training ------#
img.array = array(dim = c(nnn, ncr, ncr))
img.list = list()
for(i in 1:nrow(train_dat.id)){
  img.array[i,,] = matrix(train_dat.id$Z1[i,], ncr, ncr)
  img.list[[i]] = matrix(train_dat.id$Z1[i,], ncr, ncr)
}

if(method == 'FPLS'){
  sup_fun = FPLS_img(type= "fourier", img.list= img.list, npc = npc, ncr = ncr)
  null_obj = coxph(Surv(time,event) ~ 1, data = tdat.id, x = T, y = T) 
  
  score_names = c()
  for(q in 1:npc){
    tname = paste("score", as.character(q), sep = "")
    score_names = c(score_names, tname)
  }
  fmla = as.formula(paste("Surv(time,event) ~ ", paste(score_names, collapse= "+")))
  fitC = survfit(Surv(time, (1-event))~1, tdat.id,  type=c("kaplan-meier"))
  KM_c = summary(fitC,times=tdat.id$time,extend=TRUE) #KM_c is now ordered by time
  null.res = NULL 
  for(i in 1:nrow(tdat.id)){
    if(tdat.id[i,]$event == 1){
      idx = which(KM_c$time == tdat.id[i,]$time)
      interim = log(tdat.id[i,]$time)/KM_c$surv[idx] 
    }else{interim = 0}
    null.res = c(null.res, interim)
  }
  lam = lam_picked
  x = sup_fun[[2]]
  y.init = null.res 
  f.obj = fdata(x)
  desMat = sup_fun[[1]]
  plsCOX  = fdata2pls(fdataobj = f.obj, y = y.init, norm = F, ncomp = npc, lambda = lam)
  score_sup = x %*% t(plsCOX$rotation$data)
  
}else{
  sup_fun = sfpca_img(type = "fourier", observed = img.list, train_dat.id = train_dat.id,
                      theta = tune_picked, npc = npc, ncr = ncr)
  
  o_tfs = sup_fun[[1]]
  
  #score
  score_sup = NULL
  for(i in 1:nrow(train_dat.id)){
    score_sup = rbind(score_sup, as.vector(t(o_tfs) %*% as.vector(img.list[[i]])))
  }
  
}

score_names = c()
for(q in 1:(ncol(score_sup))){
  tname = paste("score", as.character(q), sep = "")
  score_names = c(score_names, tname)
}
tdat.id = tdat.id[,c(1:4)]
tdat.id = cbind(tdat.id, score_sup)

colnames(tdat.id)[(ncol(tdat.id) - ncol(score_sup) + 1) : ncol(tdat.id)] = score_names

fmla = as.formula(paste("Surv(time,event) ~ ", paste(score_names, collapse= "+")))

fitted_obj = coxph(fmla, data = tdat.id, x = T, y = T) 

#------------------------- TEST ---------------------------------------

img.list.test = list()
for(i in 1:nrow(test_dat.id)){
  img.list.test[[i]] = matrix(test_dat.id$Z1[i,], ncr, ncr)
}


if(method == 'FPLS'){
  XB_test = NULL
  for(i in 1:nrow(test_dat.id)){
    XB_test = rbind(XB_test, as.vector(t(desMat) %*% as.vector(img.list.test[[i]])))
  }
  score_sup_test = XB_test %*% t(plsCOX$rotation$data)
  
}else{
  score_sup_test = NULL
  for(i in 1:nrow(test_dat.id)){
    score_sup_test = rbind(score_sup_test, as.vector(t(o_tfs) %*% as.vector(img.list.test[[i]])))
  }
}


test_dat.id = test_dat.id[,c(1:4)]

score_names = c()
for(q in 1:(ncol( score_sup_test))){
  tname = paste("score", as.character(q), sep = "")
  score_names = c(score_names, tname)
}

test_dat.id = cbind(test_dat.id, score_sup_test)

colnames(test_dat.id)[(ncol(test_dat.id) - ncol(score_sup_test) + 1) : ncol(test_dat.id)] = score_names

train.prob = predict(fitted_obj)
surv.prob = predict(fitted_obj, newdata = test_dat.id)

Surv.train = Surv(tdat.id$time, tdat.id$event)
Surv.test = Surv(test_dat.id$time, test_dat.id$event)

times = seq(0, 15, by = 0.1)

# cumulative AUC 
auc.sFPCA = AUC.uno(Surv.train, Surv.test, surv.prob, times)$iauc
print(auc.sFPCA)



