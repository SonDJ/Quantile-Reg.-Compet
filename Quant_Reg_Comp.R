library(cmprskQR) ; library(survival)

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
  if(dist=='normal') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qnorm(unif, 0, 1)), exp(-rep(1,length(status))+z1-z2+qnorm(unif, 0, 1)))
  else if(dist=='logistic') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qlogis(unif, 0, 1)), exp(-rep(1,length(status))+z1-z2+qlogis(unif, 0, 1)))
  else if(dist=='cauchy') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qcauchy(unif, 0, 1)), exp(-rep(1,length(status))+z1-z2+qcauchy(unif, 0, 1)))
  else if(dist=='weibull') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qweibull(unif, 0.5, 1)), exp(-rep(1,length(status))+z1-z2+qweibull(unif, 0.5, 1)))
  }
  return(time_vec)
}

gamma=function(tau, n, obs, status, covariate, beta){
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)$surv
  if(km.cens[n]==0) km.cens[n]=km.cens[n-1]
  sum.mat.1=list()
  for(i in 1:n){
    sum.mat.1[[i]]=outer(covariate[i,], covariate[i,])*as.numeric((ifelse(log(obs[i])-(covariate%*%beta)[i]<=0 & status[i]==1, 1, 0)/km.cens[i]-tau))^2
  }
  
  l=list() ; ll=list() ; divide=c()
  for(i in 1:n){
    divide[i]=sum(ifelse(obs[i]<=obs, 1, 0))
    sub.list=list()
    for(j in 1:n){
      sub.list[[j]]=covariate[j,]*as.numeric(ifelse(obs[j]>=obs[i], 1, 0)*ifelse(obs[j]<=exp((covariate%*%beta)[j]) & status[j]==1, 1, 0)/(km.cens[j]*divide[i]))
    }
    ll[[i]]=Reduce('+', sub.list)
    l[[i]]=ifelse(status[i]==0, 1, 0)*outer(ll[[i]], ll[[i]])
  }
  return(1/n*Reduce('+', sum.mat.1)-1/n*Reduce('+', l))
}

A=function(n, obs, status, covariate, beta, sigma){
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)$surv
  if(km.cens[n]==0) km.cens[n]=km.cens[n-1]
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=as.numeric(ifelse(status[i]==1, 1, 0)/km.cens[i]*dnorm(-(log(obs[i])-(covariate%*%beta)[i])/sqrt(as.numeric(t(covariate[i,])%*%sigma%*%covariate[i,])), 0, 1))/as.numeric(sqrt(t(covariate[i,])%*%sigma%*%covariate[i,]))*outer(covariate[i,], covariate[i,])
  }
  return(1/n*Reduce('+', sum.list))
}

smooth.est.eq=function(tau, n, obs, status, covariate, beta, sigma){
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)$surv
  if(km.cens[n]==0) km.cens[n]=km.cens[n-1]
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=covariate[i,]*as.numeric(ifelse(status[i]==1, 1, 0)/km.cens[i]*pnorm(-(log(obs[i])-(covariate%*%beta)[i])/sqrt(as.numeric(t(covariate[i,])%*%sigma%*%covariate[i,])), 0, 1)-tau)
  }
  return(1/n*Reduce('+', sum.list))
}

NewtonRaphson=function(taus, ns, obss, statuss, covariates, betas){
  sols=as.matrix(betas) ; cov=1/ns*diag(rep(1, ncol(covariates)), ncol(covariates)) ; i=1
  repeat{
    new.beta=sols[,i]-solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i], sigma = cov))%*%as.matrix(smooth.est.eq(tau = taus, n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i], sigma = cov))
    sols=cbind(sols, new.beta)
    new.sigma=1/ns*solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i], sigma = cov))%*%gamma(tau = taus, n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i])%*%solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i], sigma = cov))
    cov=new.sigma
    i=i+1
    if(all(sapply(abs(sols[,i]-sols[,i-1]), FUN=function(j) {I(j<10^-5)}))) break
  }
  return(list(sols[,ncol(sols)], cov))
}

#m-repetitions of N-sample simulation with given Tau and specific distribution

simulation_function=function(m, N, Dist, L, Tau){
  reg.est=list() ; cov.est=list()
  for(i in 1:m){
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(time>cens, cens, time)
    Eps=ifelse(Obs==cens, 0, Eps)
    pfcmp=crrQR(ftime = log(Obs), fstatus = Eps, X = model.matrix(~Z1+Z2)[,-1], tau.range = c(Tau, Tau))
    a=NewtonRaphson(taus = Tau, ns = N, obss = Obs, statuss = Eps, covariates = cbind(rep(1,N), Z1, Z2), betas = as.vector(pfcmp$beta.seq))
    reg.est[[i]]=a[[1]] ; cov.est[[i]]=a[[2]]
  }
  return(list(1/m*Reduce('+', reg.est), 1/m*Reduce('+', cov.est)))
}

normal.model.300=simulation_function(m = 100, N = 300, Dist = 'normal', L = 4.6, Tau = 0.2)

#no consideration of covariance matrix

Naive.NewtonRaphson=function(taus, ns, obss, statuss, covariates, betas){
  sols=as.matrix(betas) ; cov=1/ns*diag(rep(1, ncol(covariates)), ncol(covariates)) ; i=1
  repeat{
    new.beta=sols[,i]-solve(A(n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i], sigma = cov))%*%as.matrix(smooth.est.eq(tau = taus, n = ns, obs = obss, status = statuss, covariate = covariates, beta = sols[,i], sigma = cov))
    sols=cbind(sols, new.beta)
    i=i+1
    if(all(sapply(abs(sols[,i]-sols[,i-1]), FUN=function(j) {I(j<1/(10^5))}))) break
  }
  return(list(sols[,ncol(sols)], cov))
}

naive.n.300=Naive.NewtonRaphson(taus = Tau, ns = N, obss = obs, statuss = eps, covariates = cbind(rep(1,N), Z1, Z2), betas = as.vector(pfcmp$beta.seq))

#Using Equation Solver

beta.smooth.est.eq=function(beta){
  tau=0.2; n=300; obs=n.obs.300; status=eps.300; covariate=cbind(rep(1,300),z1.300,z2.300); sigma=diag(rep(1/300,3),3)
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)$surv
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=covariate[i,]*as.numeric(ifelse(status[i]==1, 1, 0)/km.cens[i]*pnorm(-(log(obs[i])-(covariate%*%beta)[i])/sqrt(t(covariate[i,])%*%sigma%*%covariate[i,]), 0, 1)-tau)
  }
  return(1/n*Reduce('+', sum.list))
}

library(rootSolve)

m=multiroot(f = beta.smooth.est.eq, start = c(-1.4,1,0.75))
