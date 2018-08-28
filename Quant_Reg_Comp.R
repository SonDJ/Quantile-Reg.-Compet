library(survival) ; library(MASS) ; library(copula)
na.omit.list=function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
inverse.check=function(m) class(try(solve(m),silent=T))=="matrix"
quadform=function(A,x) {return(colSums(x * (A %*% x)))}
norm_vec=function(x) return(sqrt(sum(x^2)))
quadform.vec=function(Z, sig) {return(diag(Z%*%sig%*%t(Z)))}
outer_func=function(x) {return(outer(X = x, Y = x))}

censoring_rate=function(N, Dist, L, numrep=1e+03){
  num.cens=c()
  for(i in 1:numrep){
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2, p0 = 0.7, p1 = 0.8)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    num.cens[i]=length(which(Eps==0))
  }
  return(mean(num.cens)/N)
}

censoring_rate_new=function(N, L, numrep=1e+03){
  num.cens=c()
  for(i in 1:numrep){
    Z1=runif(n = N, min = 0, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2, p0 = 0.9, p1 = 0.7)
    time=exp(sampling_func_new(status = Eps.fail, z1 = Z1, z2 = Z2))
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    num.cens[i]=length(which(Eps==0))
  }
  return(mean(num.cens)/N)
}

cause_sampling_func=function(z,p0,p1){
  sample_vec=c()
  for(i in 1:length(z)){
    if(z[i]==0){
      sample_vec[i]=sample(x = c(1,2), size = 1, replace = T, prob = c(p0, 1-p0))
    }
    else{
      sample_vec[i]=sample(x = c(1,2), size = 1, replace = T, prob = c(p1, 1-p1))
    }
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

sampling_func_new=function(status, z1, z2){
  time_vec=c()
  for(i in 1:length(status)){
    if(status[i]==1){
      time_vec[i]=mvrnorm(n = 1, mu = c(0.5*z1[i]-0.5*z2[i], 0.5*z2[i]), Sigma = matrix(c(1,0.9,0.9,1), 2))[1]
    }
    else{
      time_vec[i]=mvrnorm(n = 1, mu = c(0.5*z1[i]-0.5*z2[i], 0.5*z2[i]), Sigma = matrix(c(1,0.9,0.9,1), 2))[2]
    }
  }
  return(time_vec)
}

A=function(obs, status, covariate, beta, sigma, case.sample=NULL, cohort.sample=NULL, other.weight=NULL){
  n=length(status)
  if(is.null(case.sample)==T) {
    q=1
    case.weight=rep(0, n)
  } else {
    q=length(status[which(case.sample==1)])/(length(which(status==1))-length(which(status[as.logical(cohort.sample)]==1)))
    case.weight=ifelse(status==1, 1, 0)*(1-cohort.sample)*(case.sample/q-1)
  }
  
  if(is.null(cohort.sample)==T){
    p=1
    cc.weight=rep(1, n)
  } else {
    p=length(which(as.logical(cohort.sample)==T))/n
    cc.weight=ifelse(status==1, 1, 0)+ifelse(status==1, 0, 1)*cohort.sample/p
  }
  
  if(is.null(other.weight)){
    other.weight=rep(0, n)
  } else {
    other.weight=other.weight-rep(1,n)
  }
  
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) 
  obs=obs[obs.per] 
  case.weight=case.weight[obs.per] 
  cc.weight=cc.weight[obs.per]
  status=status[obs.per] 
  covariate=covariate[obs.per,] 
  other.weight=other.weight[obs.per]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  
  dnorm.vec=dnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)/sqrt(quadvec)
  weight=ifelse(status==1, 1/km.cens, 0)
  
  return(1/n*t(covariate)%*%diag((cc.weight+case.weight+other.weight)*weight*dnorm.vec)%*%covariate)
}

smooth.gamma=function(beta, tau, obs, status, covariate, sigma, cohort.sample=NULL, case.sample=NULL){
  n=length(obs)
  if(is.null(case.sample)==T) {
    q=1
    case.weight=rep(0, n)
  } else {
    q=length(which(case.sample==1))/(length(which(status==1))-length(which(status[as.logical(cohort.sample)]==1)))
    case.weight=ifelse(status==1, 1, 0)*(1-cohort.sample)*(case.sample/q-1)
  }
  
  if(is.null(cohort.sample)==T){
    p=1
    cc.weight=rep(1, n)
  } else {
    p=length(which(cohort.sample==1))/n
    cc.weight=ifelse(status==1, 1, 0)+ifelse(status==1, 0, 1)*cohort.sample/p
  }
  
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) 
  obs=obs[obs.per] 
  case.weight=case.weight[obs.per]
  cc.weight=cc.weight[obs.per]
  status=status[obs.per] 
  covariate=covariate[obs.per,]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  
  mat=sweep(covariate, MARGIN = 1, (weight.vec*pnorm.vec-rep(tau, length(status))), '*')
  
  cen.vec=ifelse(status==0, 1, 0)
  l=matrix(rep(0, length(status)*ncol(covariate)), nrow = length(status))
  for(i in 1:length(status)){
    l[i,]=colSums(sweep(covariate, MARGIN = 1, ifelse(obs>=obs[i], 1, 0)*pnorm.vec*weight.vec, FUN = '*'))/sum(ifelse(obs>=obs[i], 1, 0))*ifelse(status[i]==0, 1, 0)
  }
  
  return(1/n*t(mat+l)%*%diag(cc.weight+case.weight)%*%(mat+l)+(1/n)*((1-p)/p)*t(mat)%*%diag(ifelse(status==1, 0, 1)*cc.weight)%*%mat+(1/n)*((1-q)/q)*(1-p)*t(mat+l)%*%diag(ifelse(status==1, 1, 0)*case.sample/q)%*%(mat+l))
}

smooth.gamma.stratified=function(beta, tau, obs, status, covariate, sigma, strata_list=NULL){
  n=length(obs)
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs)
  obs=obs[obs.per] 
  status=status[obs.per] 
  covariate=covariate[obs.per,]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  
  mat=sweep(covariate, MARGIN = 1, (weight.vec*pnorm.vec-rep(tau, length(status))), '*')
  
  cen.vec=ifelse(status==0, 1, 0)
  l=matrix(rep(0, length(status)*ncol(covariate)), nrow = length(status))
  for(i in 1:length(status)){
    l[i,]=colSums(sweep(covariate, MARGIN = 1, ifelse(obs>=obs[i], 1, 0)*pnorm.vec*weight.vec, FUN = '*'))/sum(ifelse(obs>=obs[i], 1, 0))*ifelse(status[i]==0, 1, 0)
  }
  
  if(is.null(strata_list)){
    return(1/n*t(mat+l)%*%(mat+l))
  } else {
    str_mat=strata_list[[1]][obs.per, ] ; str_weight=strata_list[[2]][obs.per]
    str_num_vec=as.vector(colSums(str_mat==1))
    alpha=str_num_vec/n
    p_vec=str_num_vec/(strata_list[[3]])-1
    strata_ind_list=lapply(apply(str_mat, 2, list), unlist)
    V_list=lapply(strata_ind_list, function(x){
      1/n*t(mat+l)%*%diag(str_weight*x)%*%(mat+l)
    })
    V_list_prod=Map(f = "*", V_list, as.list(p_vec*alpha+alpha*(1-alpha)))
    
    return(1/n*t(mat+l)%*%diag(str_weight)%*%(mat+l)+Reduce('+', V_list_prod))
  }
}

smooth.est.eq=function(beta, tau, obs, status, covariate, sigma, cohort.sample=NULL, case.sample=NULL, other.weight=NULL){
  n=length(obs)
  if(is.null(case.sample)==T) {
    q=1
    case.weight=rep(0, n)
  } else {
    q=length(status[which(case.sample==1)])/(length(which(status==1))-length(which(status[as.logical(cohort.sample)]==1)))
    case.weight=ifelse(status==1, 1, 0)*(1-cohort.sample)*(case.sample/q-1)
  }
  
  if(is.null(cohort.sample)==T){
    p=1
    cc.weight=rep(1, n)
  } else {
    p=length(which(as.logical(cohort.sample)==T))/n
    cc.weight=ifelse(status==1, 1, 0)+ifelse(status==1, 0, 1)*cohort.sample/p
  }
  
  if(is.null(other.weight)){
    other.weight=rep(0, n)
  } else {
    other.weight=other.weight-rep(1,n)
  }
  
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs)
  obs=obs[obs.per]
  case.weight=case.weight[obs.per] 
  cc.weight=cc.weight[obs.per]
  status=status[obs.per] 
  covariate=covariate[obs.per,]
  other.weight=other.weight[obs.per]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  
  return(1/n*colSums(sweep(covariate, MARGIN = 1, (cc.weight+case.weight+other.weight)*(weight.vec*pnorm.vec-rep(tau, length(status))), '*')))
}

#MB method
boot.smooth.est.eq=function(eta, beta, tau, obs, status, covariate, sigma, cohort.sample=NULL, case.sample=NULL){
  n=length(obs)
  if(is.null(case.sample)==T) {
    q=1
    case.weight=rep(0, n)
  } else {
    q=length(status[which(case.sample==1)])/(length(which(status==1))-length(which(status[as.logical(cohort.sample)]==1)))
    case.weight=ifelse(status==1, 1, 0)*(1-cohort.sample)*(case.sample/q-1)
  }
  
  if(is.null(cohort.sample)==T){
    p=1
    cc.weight=rep(1, n)
  } else {
    p=length(which(as.logical(cohort.sample)==T))/n
    cc.weight=ifelse(status==1, 1, 0)+ifelse(status==1, 0, 1)*cohort.sample/p
  }
  
  if(is.null(other.weight)){
    other.weight=rep(0, n)
  } else {
    other.weight=other.weight-rep(1,n)
  }
  
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs)
  obs=obs[obs.per]
  case.weight=case.weight[obs.per] 
  cc.weight=cc.weight[obs.per]
  status=status[obs.per] 
  covariate=covariate[obs.per,]
  other.weight=other.weight[obs.per]
  eta=eta[obs.per]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  
  return(1/n*colSums(sweep(covariate, MARGIN = 1, eta*(cc.weight+case.weight+other.weight)*(weight.vec*pnorm.vec-rep(tau, length(status))), '*')))
}

#Case-cohort design
CCHD=function(rate, data, status, case.sampling=NULL){
  all.cases=which(status==1)
  sub.coh.ind=sample(x = 1:nrow(data), size = floor(nrow(data)*rate), replace = F)
  
  if(is.null(case.sampling)==T) {
    case.coh.ind=union(all.cases, sub.coh.ind)
    case.ind=NULL
  } else {
    case.non.cohort=setdiff(all.cases, sub.coh.ind)
    case.ind=sample(x = case.non.cohort, size = floor(length(case.non.cohort)*case.sampling), replace = F)
    case.coh.ind=union(case.ind, sub.coh.ind)
  }
  
  c.vec=ifelse(1:nrow(data) %in% case.ind, 1, 0)

  return(list(ifelse(1:nrow(data) %in% sub.coh.ind, 1, 0), c.vec))
}

stratified_sampling=function(str.var, size){
  str_class=unique(str.var)[order(unique(str.var))]
  str_sample_list=list()
  strata_matrix=matrix(rep(0, length(str.var)*length(str_class)), nrow = length(str.var))
  weight_matrix=matrix(rep(0, length(str.var)*length(str_class)), nrow = length(str.var))
    for(i in 1:length(str_class)){
    strata_matrix[,i]=as.numeric(str.var %in% str_class[i])
    str_sample_list[[i]]=sample(x = which(strata_matrix[,i]==1), size = size, replace = F)
    weight_matrix[,i]=as.numeric(1:length(str.var) %in% str_sample_list[[i]])*strata_matrix[,i]/(length(str_sample_list[[i]])/length(which(strata_matrix[,i]==1)))
    }
  colnames(strata_matrix)=str_class
  return(list(strata_matrix, rowSums(weight_matrix), sapply(str_sample_list, FUN = length)))
}

# Iterative Method
Iter_simulation_gamma=function(m, N, Dist, L, Tau, Cohort=F, Case.Sample=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    P0=0.7 ; P1=0.8
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    COVAR=cbind(rep(1, N), Z1, Z2)
    if(isFALSE(Case.Sample)) CS=NULL else CS=0.8
    CCH_desig=CCHD(rate = 0.2, data = cbind(rep(1, N), Z1, Z2), status = Eps, case.sampling = CS)
    if(Cohort==T) CCH_status=CCH_desig[[1]] else CCH_status=NULL
    if(Case.Sample==T) Case_status=CCH_desig[[2]] else Case_status=NULL
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(-1+qnorm(Tau/P0), 1, 1+qnorm(Tau/P1)-qnorm(Tau/P0))
      repeat.beta=list(in.vec) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status))==F|any(is.nan(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status)))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status), smooth.est.eq(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, cohort.sample = CCH_status, case.sample = Case_status)))
          
          V.cov=smooth.gamma(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, cohort.sample = CCH_status, case.sample = Case_status)
          
          repeat.cov[[j]]=1/N*solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status), V.cov)%*%solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status))
          if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-2)|length(repeat.beta)==50) break else j=j+1
        }
      }
      if(any(is.na(repeat.beta[[length(repeat.beta)]]))==T|norm_vec(repeat.beta[[length(repeat.beta)]])>10) {m.sol[[i]]=NULL ; boot.cov[[i]]=NULL}
      else{
        m.sol[[i]]=repeat.beta[[length(repeat.beta)]]
        boot.cov[[i]]=repeat.cov[[length(repeat.cov)]]
      }
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(m.sol)*Reduce('+', boot.cov)
  return(list(1/length(m.sol)*Reduce('+', m.sol), sample.cov.mat, boot.mat))
}

Iter_simulation_gamma_new=function(m, N, L, Tau, Cohort=F, Case.Sample=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = 0, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=exp(sampling_func_new(status = Eps.fail, z1 = Z1, z2 = Z2))
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    COVAR=cbind(rep(1, N), Z1, Z2)
    if(isFALSE(Case.Sample)) CS=NULL else CS=0.5
    CCH_desig=CCHD(rate = 0.2, data = cbind(rep(1, N), Z1, Z2), status = Eps, case.sampling = CS)
    if(Cohort==T) CCH_status=CCH_desig[[1]] else CCH_status=NULL
    if(Case.Sample==T) Case_status=CCH_desig[[2]] else Case_status=NULL
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(qnorm(Tau/P0), 0.5, -0.5+qnorm(Tau/P1)-qnorm(Tau/P0))
      repeat.beta=list(in.vec) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status))==F|any(is.nan(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status)))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status), smooth.est.eq(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, cohort.sample = CCH_status, case.sample = Case_status)))
          
          V.cov=smooth.gamma(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, cohort.sample = CCH_status, case.sample = Case_status)
          
          repeat.cov[[j]]=1/N*solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status), V.cov)%*%solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status))
          if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-2)|length(repeat.beta)==50) break else j=j+1
        }
      }
      if(any(is.na(repeat.beta[[length(repeat.beta)]]))==T|norm_vec(repeat.beta[[length(repeat.beta)]])>10) {m.sol[[i]]=NULL ; boot.cov[[i]]=NULL}
      else{
        m.sol[[i]]=repeat.beta[[length(repeat.beta)]]
        boot.cov[[i]]=repeat.cov[[length(repeat.cov)]]
      }
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(m.sol)*Reduce('+', boot.cov)
  return(list((1/length(m.sol)*Reduce('+', m.sol)-in.vec), sample.cov.mat, boot.mat))
}

Iter_simulation_stratified=function(m, N, L, Tau, stratify=FALSE){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = 0, max = 1)
    Copula=rCopula(n = N, copula = normalCopula(0.8, dim = 2))
    Z2=qbinom(p = Copula[,1], size = 1, prob = 0.5)
    Z3=qbinom(p = Copula[,2], size = 1, prob = 0.8)
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=exp(sampling_func_new(status = Eps.fail, z1 = Z1, z2 = Z2))
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    COVAR=cbind(rep(1,N), Z1, Z2)
    if(isFALSE(stratify)==T) Str_list=NULL else Str_list=stratified_sampling(str.var = Z3, size = length(which(Z3==0))*0.8)
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(qnorm(Tau/P0), 0.5, -0.5+qnorm(Tau/P1)-qnorm(Tau/P0))
      repeat.beta=list(in.vec) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]))==F|any(is.nan(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]])))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]), smooth.est.eq(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, other.weight = Str_list[[2]])))
          
          V.cov=smooth.gamma.stratified(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, strata_list = Str_list)
          
          repeat.cov[[j]]=1/N*solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]), V.cov)%*%solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]))
          if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-2)|length(repeat.beta)==50) break else j=j+1
        }
      }
      if(any(is.na(repeat.beta[[length(repeat.beta)]]))==T|norm_vec(repeat.beta[[length(repeat.beta)]])>10) {m.sol[[i]]=NULL ; boot.cov[[i]]=NULL}
      else{
        m.sol[[i]]=repeat.beta[[length(repeat.beta)]]
        boot.cov[[i]]=repeat.cov[[length(repeat.cov)]]
      }
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(m.sol)*Reduce('+', boot.cov)
  return(list((1/length(m.sol)*Reduce('+', m.sol)-in.vec), sample.cov.mat, boot.mat))
}

Iter_simulation=function(m, B, N, L, Tau, Cohort=F, Case.Sample=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = 0, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=exp(sampling_func_new(status = Eps.fail, z1 = Z1, z2 = Z2))
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    if(isFALSE(Case.Sample)) CS=NULL else CS=0.5
    CCH_desig=CCHD(rate = 0.2, data = cbind(rep(1, N), Z1, Z2), status = Eps, case.sampling = CS)
    if(Cohort==T) COVAR=CCH_desig[[1]] else COVAR=cbind(rep(1, N), Z1, Z2)
    if(Cohort==T) CCH_status=CCH_desig[[2]] else CCH_status=NULL
    Case_status=CCH_desig[[3]]
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(-1+qnorm(Tau/P0), 1, 1+qnorm(Tau/P1)-qnorm(Tau/P0))
      repeat.beta=list(in.vec) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status))==F|any(is.nan(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status)))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status), smooth.est.eq(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, cohort.sample = CCH_status, case.sample = Case_status)))
          
          boot.list=list()
          for(k in 1:B){
            Eta=rexp(n = N, rate = 1)
            boot.list[[k]]=boot.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], eta = Eta, cohort.sample = CCH_status, case.sample = Case_status)
          }
          
          boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
          boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
          
          V.cov=1/(B-1)*Reduce('+', boot_matrix_sum)
          
          repeat.cov[[j]]=solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status), V.cov)%*%solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], cohort.sample = CCH_status, case.sample = Case_status))
          if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-2)|length(repeat.beta)==50) break else j=j+1
        }
      }
      if(any(is.na(repeat.beta[[length(repeat.beta)]]))==T|norm_vec(repeat.beta[[length(repeat.beta)]])>10) {m.sol[[i]]=NULL ; boot.cov[[i]]=NULL}
      else{
        m.sol[[i]]=repeat.beta[[length(repeat.beta)]]
        boot.cov[[i]]=repeat.cov[[length(repeat.cov)]]
      }
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(m.sol)*Reduce('+', boot.cov)
  return(list(1/length(m.sol)*Reduce('+', m.sol), sample.cov.mat, boot.mat))
}
