sfpca_img = function(type, observed,train_dat.id,theta,lambda,npc,timepts=NULL,ncr){
  
  N = nrow(train_dat.id)
  
  stopifnot(all(c("MFPCA","dplyr", "parallel")%in%rownames(installed.packages()))) 
  
  surv_dat = train_dat.id
  surv_log = surv_dat
  
  surv_log$logTime = log(surv_dat$time) 
  fitC = survfit(Surv(logTime, (1-event))~1, surv_log,  type=c("kaplan-meier"))
  KM_c = summary(fitC,times=surv_log$logTime,extend=TRUE) #KM_c is now ordered by time
  surv_dat = surv_log[order(surv_log$logTime),]
  surv_dat$ipcw = surv_dat$event/KM_c$surv
  surv_dat$ipcw[is.nan(surv_dat$ipcw)] = 0
  Y_bar = mean(surv_dat$ipcw * surv_dat$logTime, na.rm = T)
  surv_dat$Y_de_mean = surv_dat$logTime - Y_bar
  surv_dat = surv_dat[match(surv_log$id, surv_dat$id),]
  train_y = surv_dat$Y_de_mean * surv_dat$ipcw 
  
  if (is.null(timepts)) {
    timepts1=1:ncr
    timepts2=1:ncr
  }
  
  if(type == "bspline"){
    norder = 4
    basis1=create.bspline.basis(rangeval=c(1,ncr), norder=norder)
    basis2=create.bspline.basis(rangeval=c(1,ncr), norder=norder)
    
  }
  
  if(type == 'fourier'){
    basis1=create.fourier.basis(rangeval=c(1,ncr))
    basis2=create.fourier.basis(rangeval=c(1,ncr))
  }
  
  coord <- expand.grid(x = timepts1, y = timepts2)
  
  globalxmat <<- lapply(1:nrow(coord),function(t){
    
    evalb1 = eval.basis(coord[t,2], basis1) 
    evalb2 = eval.basis(coord[t,1], basis2)
    c(t(evalb1) %*% evalb2) 
  }
  )%>%do.call(rbind,.)
  
  dim(globalxmat)
  
  desMat = globalxmat
  
  system.time(coef.save <- mclapply(1:N, function(f) {
    lm.fit(y = as.vector(observed[[f]]), x =  desMat)$coefficients 
  }, mc.cores = 6))
  
  scores <- (do.call(rbind, coef.save))
  
  B <- aperm(array(desMat, c(ncr, ncr, 
                             ncol(scores))), c(3, 1, 2)) 
  
  mya = array(dim = c(length(observed), ncr, ncr))
  for(i in 1){
    mya[i,,] = matrix(observed[[i]], ncr, ncr)
  }
  
  g <- funData(list(c(1:ncr), c(1:ncr)), mya) 
  
  W = MFPCA:::calcBasisIntegrals(B, 2, 
                                 g@argvals)
  
  S = t(scores)
  dim(W)
  dim(S)
  
  maty = matrix(rep(train_y,each=nrow(S)),nrow=nrow(S))
  M = rowSums(maty*W%*%S)
  MM = as.matrix(M)%*%t(as.matrix(M))
  sqrM = function (X) 
  {
    EX <- eigen(X)
    VX <- EX$values
    QX <- EX$vectors
    YX <- QX %*% diag(1/sqrt(VX)) %*% t(QX)
    return(YX)
  }
  
  U =theta/length(train_y)*W%*%S%*%t(S)%*%W+(1-theta)*MM/(length(train_y)^2) 
  
  G = W 
  
  halfG_inv  = sqrM(G)
  
  tM = t(halfG_inv)%*%U%*%halfG_inv
  eigen_res = eigen(tM) 
  
  fd_list = lapply(1:npc, function(ipc){
    #cat(ipc)
    coef_pc= halfG_inv%*%as.matrix(eigen_res$vectors[,ipc]) 
    
  })
  
  basis = globalxmat  
  sup_basis = NULL
  for(k in 1:npc){
    kth = matrix(fd_list[[k]], nrow = 1) %*% t(basis) 
    sup_basis = cbind(sup_basis, t(kth))
  }
  
  return(list(sup_basis, eigen_res$values[1:npc]))
}


cond.prob = function(model, newdata, Tstart, Tpred){
  risk.Tstart = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tstart)$surv)
  risk.Tpred = as.numeric(summary(survfit(model, newdata = newdata, se.fit = F, conf.int = F), times = Tpred)$surv)
  
  return(risk.Tpred/risk.Tstart)
}

