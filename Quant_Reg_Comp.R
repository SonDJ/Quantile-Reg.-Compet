library(cmprskQR) ; library(survival) ; library(rlist) ; library(quantreg)

set.seed(100)
z1.300=runif(n = 300, min = -1, max = 1) ; z1.500=runif(n = 500, min = -1, max = 1)
z2.300=rbinom(n = 300, size = 1, prob = 0.5) ; z2.500=rbinom(n = 500, size = 1, prob = 0.5)

cause_sampling_func=function(z){
  sample_vec=c()
  for(i in 1:length(z)){
    sample_vec[i]=sample(x = c(1,2), size = 1, replace = T, prob = c(ifelse(z[i]==1, 0.8, 0.6), 1-ifelse(z[i]==1, 0.8, 0.6)))
  }
  return(sample_vec)
}

set.seed(1000)
eps.300=cause_sampling_func(z2.300) ; eps.500=cause_sampling_func(z2.500)

sampling_func=function(dist, status, z1, z2){
  unif=runif(n = length(status), min = 0, max = 1)
  time_vec=c()
  
  for(i in 1:length(status)){
  if(dist=='normal') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qnorm(unif, 0, 1)), exp(-rep(1,length(status))+z1-z2+qnorm(unif, 0, 1)))
  else if(dist=='logistic') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qlogis(unif, 0, 1)), exp(-rep(1,length(status))+z1-z2+qlogis(unif, 0, 1)))
  else if(dist=='cauchy') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qcauchy(unif, 0, 1)), exp(-rep(1,length(status))+z1-z2+qcauchy(unif, 0, 1)))
  else if(dist=='weibull') time_vec=ifelse(status==1, exp(-rep(1,length(status))+z1+z2+qweibullm(unif, 0.5, 1)), exp(-rep(1,length(status))+z1-z2+qweibull(unif, 0.5, 1)))
  }
  return(time_vec)
}

gamma=function(tau, n, obs, status, covariate, beta){
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)
  sum.mat.1=list()
  for(i in 1:n){
    sum.mat.1[[i]]=outer(covariate[i,], covariate[i,])*as.numeric((ifelse(log(obs[i])-(covariate%*%beta)[i]<=0 & status[i]==1, 1, 0)/km.cens$surv[i]-tau))^2
  }
  
  l=list() ; ll=list() ; divide=c()
  for(i in 1:n){
    divide[i]=sum(ifelse(obs[i]<=obs, 1, 0))
    sub.list=list()
    for(j in 1:n){
      sub.list[[j]]=covariate[j,]*as.numeric(ifelse(obs[j]>=obs[i], 1, 0)*ifelse(obs[j]<=exp((covariate%*%beta)[j]) & status[j]==1, 1, 0)/(km.cens$surv[j]*divide[i]))
    }
    ll[[i]]=Reduce('+', sub.list)
    l[[i]]=ifelse(status[i]==0, 1, 0)*outer(ll[[i]], ll[[i]])
  }
  return(1/n*Reduce('+', sum.mat.1)-1/n*Reduce('+', l))
}

A=function(n, obs, status, covariate, beta, sigma){
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)$surv
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=as.numeric(ifelse(status[i]==1, 1, 0)/km.cens[i]*dnorm(-(log(obs[i])-(covariate%*%beta)[i])/sqrt(t(covariate[i,])%*%sigma%*%covariate[i,]), 0, 1))/as.numeric(sqrt(t(covariate[i,])%*%sigma%*%covariate[i,]))*outer(covariate[i,], covariate[i,])
  }
  return(1/n*Reduce('+', sum.list))
}

smooth.est.eq=function(tau, n, obs, status, covariate, beta, sigma){
  km.cens=survfit(Surv(time = obs, event = ifelse(status==0,1,0))~1)$surv
  sum.list=list()
  for(i in 1:n){
    sum.list[[i]]=covariate[i,]*as.numeric(ifelse(status[i]==1, 1, 0)/km.cens[i]*pnorm(-(log(obs[i])-(covariate%*%beta)[i])/sqrt(t(covariate[i,])%*%sigma%*%covariate[i,]), 0, 1)-tau)
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
    if(all(sapply(abs(sols[,i]-sols[,i-1]), FUN=function(j) {I(j<1/(10^5))}))) break
  }
  return(list(sols[,ncol(sols)], cov))
}

#n=300 cases-Normal

set.seed(300)
n.time.300=sampling_func(dist = 'normal', status = eps.300, z1 = z1.300, z2 = z2.300)

set.seed(400)
cens.300=runif(300, 0, 4.6) ; cens.500=runif(500, 0, 4.6)

n.obs.300=ifelse(n.time.300>cens.300, cens.300, n.time.300)

eps.300=ifelse(n.obs.300==cens.300, 0, eps.300)

pfcmp.300=crrQR(ftime = n.obs.300, fstatus = eps.300, X = model.matrix(~z1.300+z2.300)[,-1], tau.step = 0.2, tau.range = c(0.2, 0.4), failcode = 1, cencode = 0)

a=NewtonRaphson(taus = 0.2, ns = 300, obss = n.obs.300, statuss = eps.300, covariates = cbind(rep(1,300), z1.300, z2.300), betas = pfcmp.300$beta.seq[1,])

j=1
sols=as.matrix(pfcmp.300$beta.seq[1,])
cov=list(diag(rep(1/300, 3), 3))
repeat{
  new.beta=sols[,j]-solve(A(n = 300, obs = n.obs.300, status = eps.300, covariate = cbind(rep(1,300),z1.300,z2.300), beta = sols[,j], sigma = cov[[j]]))%*%as.matrix(smooth.est.eq(tau = 0.2, n = 300, obs = n.obs.300, status = eps.300, covariate = cbind(rep(1,300),z1.300,z2.300),beta = sols[,j], sigma = cov[[j]]))
  sols=cbind(sols, new.beta)
  new.sigma=1/300*solve(A(n = 300, obs = n.obs.300, status = eps.300, covariate = cbind(rep(1,300), z1.300, z2.300), beta = sols[,j], sigma = cov[[j]]))%*%gamma(tau = 0.2, n = 300, obs = n.obs.300, status = eps.300, covariate = cbind(rep(1,300), z1.300, z2.300), beta = sols[,j])%*%solve(A(n = 300, obs = n.obs.300, status = eps.300, covariate = cbind(rep(1,300), z1.300, z2.300), beta = sols[,j], sigma = cov[[j]]))
  cov=list.append(cov, new.sigma)
  j=j+1
  if(all(sapply(abs(sols[,j]-sols[,j-1]), FUN=function(i) {I(i<1/(10^5))}))) break
}

tau=0.2 ; n=300 ; obs=n.obs.300 ; status=eps.300 ; covariate=cbind(rep(1, 300), z1.300, z2.300) ; beta=c(-1.4, 1, 0.756)