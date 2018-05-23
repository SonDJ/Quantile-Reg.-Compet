library(cmprskQR) ; library(survival) ; library(nleqslv) ; library(MASS) ; library(rootSolve)
na.omit.list=function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

cause_sampling_func=function(z){
  sample_vec=c()
  for(i in 1:length(z)){
    sample_vec[i]=sample(x = c(1,2), size = 1, replace = T, prob = c(ifelse(z[i]==1, 0.8, 0.6), ifelse(z[i]==1, 0.2, 0.4)))
  }
  return(sample_vec)
}

sampling_func=function(dist, status, z1, z2){
  unif=runif(n = length(status), min = 0, max = 1)
  time_vec=c()
  
  for(i in 1:length(status)){
    if(dist=='normal') time_vec[i]=ifelse(status[i]==1, exp(-1+z1[i]+z2[i]+qnorm(unif[i], 0, 1)), exp(-1+z1[i]-z2[i]+qnorm(unif[i], 0, 1)))
    else if(dist=='logistic') time_vec[i]=ifelse(status[i]==1, exp(-1+z1[i]+z2[i]+qlogis(unif[i], 0, 1)), exp(-1+z1[i]-z2[i]+qlogis(unif[i], 0, 1)))
    else if(dist=='cauchy') time_vec[i]=ifelse(status[i]==1, exp(-1+z1[i]+z2[i]+qcauchy(unif[i], 0, 1)), exp(-1+z1[i]-z2[i]+qcauchy(unif[i], 0, 1)))
    else if(dist=='weibull') time_vec[i]=ifelse(status[i]==1, exp(-1+z1[i]+z2[i]+qweibull(unif[i], 0.5, 1)), exp(-1+z1[i]-z2[i]+qweibull(unif[i], 0.5, 1)))
  }
  return(time_vec)
}

A=function(n, obs, status, covariate, beta, sigma, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=c()
  for(i in 1:n){
    obs.per[i]=which(obs==SF$time[i])
  }
  if(is.null(CCH)) pn=1 else pn=CCH
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=ifelse(status[obs.per[i]]==1, 1, ifelse(any(is.na(covariate[obs.per[i],]))==F, 1/pn, 0))*ifelse(status[which(obs==SF$time[i])]==1, 1/km.cens[i], 0)*dnorm(-(log(obs[obs.per[i]])-as.numeric((covariate%*%beta)[obs.per[i]]))/sqrt(as.numeric(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],])), 0, 1)/as.numeric(sqrt(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],]))*outer(covariate[obs.per[i],], covariate[obs.per[i],])
  }
  sum.list=na.omit.list(sum.list)
  return(1/n*Reduce('+', sum.list))
}

smooth.est.eq=function(beta, tau, n, obs, status, covariate, sigma, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=match(x = SF$time, table = obs)
  if(is.null(CCH)) pn=1 else pn=CCH
  sum.list=list()
  for(i in 1:n){
    if(status[obs.per[i]]==1){
      sum.list[[i]]=as.vector(covariate[obs.per[i],])*(1/km.cens[i]*pnorm(-(log(obs[obs.per[i]])-as.numeric((covariate%*%beta)[obs.per[i]]))/sqrt(as.numeric(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],])), 0, 1)-tau)
    }
    else{
      sum.list[[i]]=-ifelse(any(is.na(covariate[obs.per[i],]))==F, 1/pn, 0)*covariate[obs.per[i],]*tau
    }
  }
  sum.list=na.omit.list(sum.list)
  return(1/n*Reduce('+', sum.list))
}

#MB method
boot.smooth.est.eq=function(tau, n, obs, status, covariate, beta, sigma, eta, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=match(x = SF$time, table = obs)
  if(is.null(CCH)) pn=1 else pn=CCH
  sum.list=list()
  for(i in 1:n){
    if(status[obs.per[i]]==1){
      sum.list[[i]]=eta[[obs.per[i]]]*as.vector(covariate[obs.per[i],])*(1/km.cens[i]*pnorm(-(log(obs[obs.per[i]])-as.numeric((covariate%*%beta)[obs.per[i]]))/sqrt(as.numeric(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],])), 0, 1)-tau)
    }
    else{
      sum.list[[i]]=-eta[[obs.per[i]]]*ifelse(any(is.na(covariate[obs.per[i],]))==F, 1/pn, 0)*covariate[obs.per[i],]*tau
    }
  }
  sum.list=na.omit.list(sum.list)
  return(1/n*Reduce('+', sum.list))
}

#Case-cohort design
CCHD=function(rate, data, status){
  location=which(status==1)
  sub.coh.ind=sample(x = 1:nrow(data), size = floor(nrow(data)*rate), replace = F)
  case.coh.ind=union(location, sub.coh.ind)
  data[setdiff(x = 1:nrow(data), y = case.coh.ind), ]=rep(NA, ncol(data))
  return(list(data, rate, sub.coh.ind))
}

#IS-MB method
ISMB_simulation_function=function(m, B, N, Dist, L, Tau, Cohort=F){
  cov.est=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(time>cens, cens, time)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    CCH_desig=CCHD(rate = 0.1, data = cbind(rep(1, N), Z1, Z2), status = Eps)
    COVAR=if(Cohort==T) CCH_desig[[1]] else cbind(rep(1, N), Z1, Z2)
    if(Cohort==T) CCH_status=CCH_desig[[2]] else CCH_status=NULL
    
    if(length(survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)$time)==N){
      pfcmp=crrQR(ftime = log(Obs), fstatus = Eps, X = model.matrix(~Z1+Z2)[,-1], tau.range = c(Tau, Tau))
      m.sol[[i]]=nleqslv(x = as.vector(pfcmp$beta.seq), fn = smooth.est.eq, method = c("Newton"), tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, sigma = cov.est, CCH = CCH_status)$x
      
      boot.list=list() ; Eta=list()
      
      for(j in 1:B){
        Eta[[j]]=rexp(n = N, rate = 1)
        boot.list[[j]]=boot.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, beta = m.sol[[i]], sigma = cov.est, eta = Eta[[j]], CCH = CCH_status)
      }
      
      boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
      boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
      V.cov=1/(B-1)*Reduce('+', boot_matrix_sum)
      A_mat=A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = m.sol[[i]], sigma = cov.est, CCH = CCH_status)
      if(rcond(A_mat)>=1e-15) solve.A = solve(A_mat) else if(rcond(A_mat)<1e-15) solve.A=ginv(X = A_mat)
      boot.cov[[i]]=solve.A%*%V.cov%*%solve.A
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/length(m.sol)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(m.sol)*Reduce('+', boot.cov)
  return(list(1/length(m.sol)*Reduce('+', m.sol), sample.cov.mat, boot.mat))
}

set.seed(300)
normal.model.com.300=ISMB_simulation_function(m = 300, B = 500, N = 300, Dist = 'normal', L = 4.2, Tau = 0.2)

set.seed(300)
normal.model.com.2500=ISMB_simulation_function(m = 2500, B = 500, N = 300, Dist = 'normal', L = 4.2, Tau = 0.2)

set.seed(300)
normal.model.cch.100.2=ISMB_simulation_function(m = 100, B = 500, N = 1500, Dist = 'normal', L = 0.5, Tau = 0.2, Cohort = T)

set.seed(300)
normal.model.cch.1000=ISMB_simulation_function(m = 1000, B = 500, N = 1500, Dist = 'normal', L = 1, Tau = 0.2, Cohort = T)
