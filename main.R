library(MFPCA);library(tidyverse);library(MASS);require(zipfR);require(fOptions)
require(Matrix);library(parallel);library(foreach);library(doMC);library(rms);library(pec)
library(fda);library(survAUC)

source("ancillary.R")

#------- one run of the simulation under simulation scenario 1 ------#

n = 400;nnn = 300
ncr = 32
npc = 1 # number of basis to pick
tune_picked = 0.01 # tunning parameter, pre-specified in this example
train_dat.id = tdat.id = read_rds('train_dat')
test_dat.id = read_rds('test_dat')

#----- training ------#
img.array = array(dim = c(nnn, ncr, ncr))
img.list = list()
for(i in 1:nrow(train_dat.id)){
  img.array[i,,] = matrix(train_dat.id$Z1[i,], ncr, ncr)
  img.list[[i]] = matrix(train_dat.id$Z1[i,], ncr, ncr)
}

sup_fun = sfpca_img(type = "fourier", observed = img.list, train_dat.id = train_dat.id,
                    theta = tune_picked, npc = npc, ncr = ncr)

o_tfs = sup_fun[[1]]

#score
score_sup = NULL
for(i in 1:nrow(train_dat.id)){
  score_sup = rbind(score_sup, as.vector(t(o_tfs) %*% as.vector(img.list[[i]])))
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

score_sup_test = NULL
for(i in 1:nrow(test_dat.id)){
  score_sup_test = rbind(score_sup_test, as.vector(t(o_tfs) %*% as.vector(img.list.test[[i]])))
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

times = seq(0.1, 15, by = 0.1)

# cumulative AUC 
auc.sFPCA = AUC.uno(Surv.train, Surv.test, surv.prob, times)$iauc
print(auc.sFPCA)



