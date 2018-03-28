library(cmprskQR) ; library(survival) ; library(nleqslv) ; library(rootSolve)

cause_sampling_func=function(z){
  sample_vec=c()
  for(i in 1:length(z)){
    sample_vec[i]=sample(x = c(1,2), size = 1, replace = T, prob = c(ifelse(z[i]==1, 0.8, 0.6), 1-ifelse(z[i]==1, 0.8, 0.6)))
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

gamma=function(tau, n, obs, status, covariate, beta){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=c()
  for(i in 1:n){
    obs.per[i]=which(obs==SF$time[i])
  }
  sum.mat.1=list()
  
  for(i in 1:n){
    sum.mat.1[[i]]=outer(covariate[obs.per[i],], covariate[obs.per[i],])*as.numeric((ifelse(log(obs[obs.per[i]])-(covariate%*%beta)[obs.per[i]]<=0 & status[which(obs==SF$time[i])]==1, 1/km.cens[i], 0)-tau))^2
  }
  
  l=list() ; ll=list() ; divide=c()
  for(i in 1:n){
    divide[i]=sum(ifelse(obs[obs.per[i]]<=obs, 1, 0))
    sub.list=list()
    for(j in 1:n){
      sub.list[[j]]=covariate[obs.per[j],]*as.numeric(ifelse(obs[obs.per[j]]>=obs[i], 1, 0)*ifelse(obs[obs.per[j]]<=exp((covariate%*%beta)[obs.per[j]]) & status[which(obs==SF$time[j])]==1, 1, 0)/(km.cens[j]*divide[i]))
    }
    ll[[i]]=Reduce('+', sub.list)
    l[[i]]=ifelse(status[obs.per[i]]==0, 1, 0)*outer(ll[[i]], ll[[i]])
  }
  return(1/n*Reduce('+', sum.mat.1)-1/n*Reduce('+', l))
}

A=function(n, obs, status, covariate, beta, sigma){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=c()
  for(i in 1:n){
    obs.per[i]=which(obs==SF$time[i])
  }
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=as.numeric(ifelse(status[which(obs==SF$time[i])]==1, 1/km.cens[i], 0)*dnorm(-(log(obs[obs.per[i]])-(covariate%*%beta)[obs.per[i]])/sqrt(as.numeric(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],])), 0, 1))/as.numeric(sqrt(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],]))*outer(covariate[obs.per[i],], covariate[obs.per[i],])
  }
  return(1/n*Reduce('+', sum.list))
}

smooth.est.eq=function(beta, tau, n, obs, status, covariate, sigma){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=c()
  for(i in 1:n){
    obs.per[i]=which(obs==SF$time[i])
  }
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=as.vector(covariate[obs.per[i],])*as.numeric(ifelse(status[which(obs==SF$time[i])]==1, 1/km.cens[i], 0)*pnorm(-(log(obs[obs.per[i]])-(covariate%*%beta)[obs.per[i]])/sqrt(as.numeric(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],])), 0, 1)-tau)
  }
  return(1/n*Reduce('+', sum.list))
}

NewtonRaphson=function(taus, ns, obss, statuss, covariates, betas){
  sols=list(betas) ; cov=1/ns*diag(rep(1, ncol(covariates)), ncol(covariates)) ; i=1
  f=function(x) class(try(solve(x),silent=T))=="matrix"
  
  repeat{
    if(f(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]), sigma = cov))==FALSE|any(is.nan(smooth.est.eq(tau = taus, n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]), sigma = cov)))) {sols[[i]]=NA ; break
    } else {
      new.beta=as.matrix(sols[[i]])-solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]), sigma = cov))%*%as.matrix(smooth.est.eq(tau = taus, n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]), sigma = cov))
      new.sigma=1/ns*solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]), sigma = cov))%*%gamma(tau = taus, n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]))%*%solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = as.matrix(sols[[i]]), sigma = cov))
      cov=new.sigma
      i=i+1
      sols[[i]]=as.vector(new.beta)
      if(all(sapply(abs(sols[[i]]-sols[[i-1]]), FUN=function(j) {I(j<10^-3)}))) break
    }
  }
  if(any(unlist(lapply(sols, is.null)))) NA else list(sols[[length(sols)]], cov)
}

#m-repetitions of N-sample simulation with given Tau and specific distribution

simulation_function=function(m, N, Dist, L, Tau){
  reg.est=list() ; cov.est=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(time>cens, cens, time)
    Eps=ifelse(Obs==cens, 0, Eps)
    pfcmp=suppressMessages(crrQR(ftime = log(Obs), fstatus = Eps, X = model.matrix(~Z1+Z2)[,-1], tau.range = c(Tau, Tau)))
    a=NewtonRaphson(taus = Tau, ns = N, obss = Obs, statuss = Eps, covariates = cbind(rep(1,N), Z1, Z2), betas = as.vector(pfcmp$beta.seq))
    reg.est[[i]]=a[[1]] ; cov.est[[i]]=a[[2]]
    i=i+1
    if(length(reg.est)==m) break
  }
  return(list(1/length(reg.est)*Reduce('+', reg.est), 1/length(reg.est)*Reduce('+', cov.est)))
}

normal.model.300=simulation_function(m = 100, N = 300, Dist = 'normal', L = 4.6, Tau = 0.2)

#MB method
MB.smooth.est.eq=function(tau, n, obs, status, covariate, beta, sigma, eta){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=c()
  for(i in 1:n){
    obs.per[i]=which(obs==SF$time[i])
  }
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=eta[obs.per[i]]*covariate[obs.per[i],]*as.numeric(ifelse(status[which(obs==SF$time[i])]==1, 1/km.cens[i], 0)*pnorm(-(log(obs[obs.per[i]])-(covariate%*%beta)[obs.per[i]])/sqrt(as.numeric(t(covariate[obs.per[i],])%*%sigma%*%covariate[obs.per[i],])), 0, 1)-tau)
  }
  return(1/n*Reduce('+', sum.list))
}

MB_simulation_function=function(m, B, N, Dist, L, Tau){
  cov.est=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list()
  for(i in 1:m){
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(time>cens, cens, time)
    Eps=ifelse(Obs==cens, 0, Eps)
    pfcmp=suppressMessages(crrQR(ftime = log(Obs), fstatus = Eps, X = model.matrix(~Z1+Z2)[,-1], tau.range = c(Tau, Tau)))
    m.sol[[i]]=nleqslv(x = as.vector(pfcmp$beta.seq), fn = smooth.est.eq, method = c("Newton"), tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), sigma = cov.est)$x
    
    boot.list=list()
    for(j in 1:B){
      Eta=rexp(n = N, rate = 1)
      boot.list[[j]]=nleqslv(x = as.vector(pfcmp$beta.seq), fn = MB.smooth.est.eq, method = c("Newton"), tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), sigma = cov.est, eta = Eta)$x
    }
    boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
    covariance_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
    boot.cov[[i]]=1/B*Reduce('+', covariance_matrix_sum)
  }
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/m*Reduce('+', covariance_matrix_sum)
  boot.mat=1/m*Reduce('+', boot.cov)
  return(list(1/length(m.sol)*Reduce('+', m.sol), sample.cov.mat, boot.mat))
}

set.seed(300)
naive.normal.model.300=MB_simulation_function(m = 100, B=100, N = 300, Dist = 'normal', L = 4.6, Tau = 0.2)

#IS-MB method
ISMB_simulation_function=function(m, B, N, Dist, L, Tau){
  cov.est=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list()
  for(i in 1:m){
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(time>cens, cens, time)
    Eps=ifelse(Obs==cens, 0, Eps)
    pfcmp=suppressMessages(crrQR(ftime = log(Obs), fstatus = Eps, X = model.matrix(~Z1+Z2)[,-1], tau.range = c(Tau, Tau)))
    m.sol[[i]]=nleqslv(x = as.vector(pfcmp$beta.seq), fn = smooth.est.eq, method = c("Newton"), tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), sigma = cov.est)$x
    
    boot.list=list()
    for(j in 1:B){
      Eta=rexp(n = N, rate = 1)
      boot.list[[j]]=MB.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), beta = as.vector(pfcmp$beta.seq), sigma = cov.est, eta = Eta)
    }
    boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
    boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
    V.cov=1/B*Reduce('+', boot_matrix_sum)
    A_mat=gradient(f = smooth.est.eq, x = as.vector(pfcmp$beta.seq), tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), sigma = cov.est)
    boot.cov[[i]]=solve(A_mat)%*%V.cov%*%solve(A_mat)
  }
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/m*Reduce('+', covariance_matrix_sum)
  boot.mat=1/m*Reduce('+', boot.cov)
  return(list(1/length(m.sol)*Reduce('+', m.sol), sample.cov.mat, boot.mat))
}

set.seed(300)
ISMB.normal.model.300=ISMB_simulation_function(m = 100, B = 100, N = 300, Dist = 'normal', L = 4.6, Tau = 0.2)