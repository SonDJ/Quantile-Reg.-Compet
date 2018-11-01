library(survival) ; library(MASS) ; library(bindata) ; library(quantreg)
na.omit.list=function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
inverse.check=function(m) class(try(solve(m),silent=T))=="matrix"
quadform=function(A,x) {return(colSums(x * (A %*% x)))}
norm_vec=function(x) return(sqrt(sum(x^2)))
quadform.vec=function(Z, sig) {return(diag(Z%*%sig%*%t(Z)))}
outer_func=function(x) {return(outer(X = x, Y = x))}

a_func=function(p,q,rho){
  return(rho*sqrt(p*q*(1-p)*(1-q))+(1-p)*(1-q))
}

corr_func=function(N, rep=1e+03){
  corr_vec=c()
  for(i in 1:rep){
    Copula=rCopula(n = N, copula = normalCopula(0.8, dim = 2))
    Z2=qbinom(p = Copula[,1], size = 1, prob = 0.5)
    Z3=qbinom(p = Copula[,2], size = 1, prob = 0.8)
    corr_vec[i]=cor(Z2,Z3)
  }
  return(mean(corr_vec))
}

censoring_rate=function(N, Dist, L, numrep=1e+03){
  num.cens=c()
  for(i in 1:numrep){
    Z1=runif(n = N, min = -1, max = 1) 
    Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2, p0 = 0.7, p1 = 0.8)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=pmin(time,cens)
    Eps=(time<=cens)*Eps.fail
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

status_prop=function(N, numrep=1e+03){
  ratio=c()
  for(i in 1:numrep){
    Z2=rbinom(n = N, size = 1, prob = 0.5)
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    ratio[i]=length(which(Eps.fail==1))/length(which(Eps.fail==2))
  }
  return(mean(ratio))
}

cause_sampling_func=function(z,p0,p1){
  n=length(z)
  pca=ifelse(z==1,p1,p0)
  ca=rbinom(n,1,prob=pca)
  cause=ifelse(ca==1,1,2)
  return(cause)
}

sampling_func=function(dist, status, z1, z2){
  n=length(status)
  beta1=c(-1, 1, 1)
  beta2=c(-1, 1, -1)
  x=cbind(rep(1,n),z1,z2)
  if(dist==1){
    T1=exp(x%*%beta1+qnorm(runif(n)))
    T2=exp(x%*%beta2+qnorm(runif(n)))
  } else if(dist==2){
    T1=exp(x%*%beta1+qlogis(runif(n)))
    T2=exp(x%*%beta2+qlogis(runif(n)))
  } else if(dist==3){
    T1=exp(x%*%beta1+qcauchy(runif(n)))
    T2=exp(x%*%beta2+qcauchy(runif(n)))
  }
  time_vec=rep(0,n)
  time_vec[status==1]=T1[status==1]
  time_vec[status==2]=T2[status==2]
  return(time_vec)
}

sampling_func_new=function(status, z1, z2){
  time_vec=c()
  for(i in 1:length(status)){
    MVNORM=mvrnorm(n = 1, mu = c(0.5*z1[i]-0.5*z2[i], 0.5*z2[i]), Sigma = matrix(c(1,0.9,0.9,1), 2))
    if(status[i]==1){
      time_vec[i]=MVNORM[1]
    }
    else{
      time_vec[i]=MVNORM[2]
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
  
  if(is.null(other.weight)==T){
    other.weight=rep(0, n)
  } else {
    other.weight=other.weight-cc.weight-case.weight
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

A.strat=function(obs, status, covariate, beta, sigma, other.weight=NULL){
  n=length(status)
  if(is.null(other.weight)==T){
    other.weight=rep(1, n)
  } else {
    other.weight=other.weight
  }
  
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) 
  obs=obs[obs.per] 
  status=status[obs.per] 
  covariate=covariate[obs.per,] 
  other.weight=other.weight[obs.per]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  
  dnorm.vec=dnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)/sqrt(quadvec)
  weight=ifelse(status==1, 1/km.cens, 0)
  
  return(1/n*t(covariate)%*%diag(other.weight*weight*dnorm.vec)%*%covariate)
}

smooth.gamma=function(beta, tau, obs, status, covariate, sigma, cohort.sample=NULL, case.sample=NULL){
  n=length(obs)
  if(is.null(case.sample)==T) {
    q=1
    case.weight=rep(0, n) ; case.sample=rep(0, n)
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
  case.weight=case.weight[obs.per] ; case.sample=case.sample[obs.per]
  cc.weight=cc.weight[obs.per] ; cohort.sample=cohort.sample[obs.per]
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
  
  l=matrix(rep(0, length(status)*ncol(covariate)), nrow = length(status))
  for(i in 1:length(status)){
    l[i,]=colSums(sweep(covariate, MARGIN = 1, ifelse(obs>=obs[i], 1, 0)*pnorm.vec*weight.vec, FUN = '*'))/sum(ifelse(obs>=obs[i], 1, 0))*ifelse(status[i]==0, 1, 0)
  }
  
  if(is.null(strata_list)==T){
    return(1/n*t(mat+l)%*%(mat+l))
  } else {
    str_mat=strata_list[[1]][obs.per, ] ; str_weight=strata_list[[2]][obs.per]
    str_num_vec=as.vector(colSums(str_mat==1))
    alpha=str_num_vec/n
    p_vec=(strata_list[[3]])/str_num_vec
    strata_ind_list=lapply(apply(str_mat, 2, list), unlist)
    V_list=list()
    for(i in 1:length(strata_ind_list)){
      V_list[[i]]=1/n*t(sweep(mat+l, MARGIN = 1, strata_ind_list[[i]], FUN = '*'))%*%diag(str_weight)%*%sweep(mat+l, MARGIN = 1, strata_ind_list[[i]], FUN = '*')-1/(n^2)*t(sweep(mat+l, MARGIN = 1, strata_ind_list[[i]]*str_weight, FUN = '*'))%*%sweep(mat+l, MARGIN = 1, strata_ind_list[[i]]*str_weight, FUN = '*')
    }
    V_list_prod=Map(f = "*", V_list, (1-p_vec)*alpha/p_vec)
    
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
    other.weight=other.weight-cc.weight-case.weight
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

smooth.est.eq.strat=function(beta, tau, obs, status, covariate, sigma, other.weight=NULL){
  n=length(obs)
  if(is.null(other.weight)){
    other.weight=rep(1, n)
  } else {
    other.weight=other.weight
  }
  
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs)
  obs=obs[obs.per]
  status=status[obs.per] 
  covariate=covariate[obs.per,]
  other.weight=other.weight[obs.per]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  
  return(1/n*colSums(sweep(covariate, MARGIN = 1, other.weight*(weight.vec*pnorm.vec-rep(tau, length(status))), '*')))
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
    str_sample_list[[i]]=sample(x = which(strata_matrix[,i]==1), size = size[i], replace = F)
    q=length(str_sample_list[[i]])/length(which(strata_matrix[,i]==1))
    weight_matrix[,i]=strata_matrix[,i]*(as.numeric(1:length(str.var) %in% str_sample_list[[i]]))/q
  }
  return(list(strata_matrix, rowSums(weight_matrix), sapply(str_sample_list, FUN = length)))
}

# Iterative Method
Iter_simulation_gamma=function(m, N, Dist, L, Tau, stratify=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; coverage=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1)
    Corr_Mat=rmvbin(n = N, margprob = c(0.5, 0.5), bincorr = matrix(c(1,0.9,0.9,1),2))
    Z2=Corr_Mat[,1]
    Z3=Corr_Mat[,2]
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=pmin(time,cens)
    Eps=(time<=cens)*Eps.fail
    New.Strata=ifelse(Eps==1, 1, ifelse(Z3==0, 2, 3))
    if(isFALSE(stratify)==T) Str_list=NULL else Str_list=stratified_sampling(str.var = New.Strata, size = c(length(which(New.Strata==1)), 150, 150))
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
      coverage[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(-1+qnorm(Tau/P0), 1, 1+qnorm(Tau/P1)-qnorm(Tau/P0))
      tryCatch({
        NR=NewtonRhapson(COVAR = cbind(Z1,Z2), Obs = Obs, Eps = Eps, Tau = Tau, In.vec = in.vec, Str_list = Str_list)
        if(norm_vec(NR[[1]])>10){
          m.sol[[i]]=NULL
          boot.cov[[i]]=NULL
          coverage[[i]]=NULL
        } else {
          m.sol[[i]]=NR[[1]]
          boot.cov[[i]]=NR[[2]]
          coverage[[i]]=as.numeric(m.sol[[i]]-qnorm(p = 0.975)*sqrt(diag(boot.cov[[i]]))<in.vec & m.sol[[i]]+qnorm(p = 0.975)*sqrt(diag(boot.cov[[i]]))>in.vec)
          
        }
      }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  coverage=Filter(Negate(is.null), coverage)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(boot.cov)*Reduce('+', boot.cov)
  return(list((1/length(m.sol)*Reduce('+', m.sol)-in.vec), sample.cov.mat, boot.mat, 1/length(coverage)*Reduce('+', coverage)))
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
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; coverage=list() ; i=1
  repeat{
    Z1=runif(n = N, min = 0, max = 1)
    Corr_Mat=rmvbin(n = N, margprob = c(0.5, 0.5), bincorr = matrix(c(1,0.9,0.9,1),2))
    Z2=Corr_Mat[,1]
    Z3=Corr_Mat[,2]
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=exp(sampling_func_new(status = Eps.fail, z1 = Z1, z2 = Z2))
    cens=runif(n = N, min = 0, max = L)
    Obs=pmin(time,cens)
    Eps=(time<=cens)*Eps.fail
    New.Strata=ifelse(Eps==1, 1, ifelse(Z3==0, 2, 3))
    COVAR=cbind(rep(1,N), Z1, Z2)
    if(isFALSE(stratify)==T) Str_list=NULL else Str_list=stratified_sampling(str.var = New.Strata, size = c(length(which(New.Strata==1)), 100, 200))
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
      coverage[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(qnorm(Tau/P0), 0.5, -0.5+qnorm(Tau/P1)-qnorm(Tau/P0))
      repeat.beta=list(in.vec) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]))==F|any(is.nan(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]])))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]), smooth.est.eq(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, other.weight = Str_list[[2]])))
          
          V.cov=smooth.gamma.stratified(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], tau = Tau, strata_list = Str_list)
          
          if(inverse.check(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]))==F|any(is.nan(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]])))==T) {repeat.cov[[j]]=NA ; break}
          
          else{
            repeat.cov[[j]]=1/N*solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]), V.cov)%*%solve(A(obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]))
          }
          
          if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-3)|length(repeat.beta)==50) break else j=j+1
        }
      }
      if(any(is.na(repeat.beta[[length(repeat.beta)]]))==T|any(is.na(repeat.cov[[length(repeat.cov)]]))==T|norm_vec(repeat.beta[[length(repeat.beta)]])>10) {m.sol[[i]]=NULL ; boot.cov[[i]]=NULL ; coverage[[i]]=NULL}
      else{
        m.sol[[i]]=repeat.beta[[length(repeat.beta)]]
        boot.cov[[i]]=repeat.cov[[length(repeat.cov)]]
        coverage[[i]]=as.numeric(m.sol[[i]]-qnorm(p = 0.975)*sqrt(diag(boot.cov[[i]]))<in.vec & m.sol[[i]]+qnorm(p = 0.975)*sqrt(diag(boot.cov[[i]]))>in.vec)
      }
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  coverage=Filter(Negate(is.null), coverage)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(boot.cov)*Reduce('+', boot.cov)
  return(list((1/length(m.sol)*Reduce('+', m.sol)-in.vec), sample.cov.mat, boot.mat, 1/length(coverage)*Reduce('+', coverage)))
}

Iter_simulation_stratified_new=function(m, N, L, Tau, stratify=FALSE){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; coverage=list() ; i=1
  repeat{
    Z1=runif(n = N, min = 0, max = 1)
    Corr_Mat=rmvbin(n = N, margprob = c(0.5, 0.5), bincorr = matrix(c(1,0.9,0.9,1),2))
    Z2=Corr_Mat[,1]
    Z3=Corr_Mat[,2]
    P0=0.9 ; P1=0.7
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=exp(sampling_func_new(status = Eps.fail, z1 = Z1, z2 = Z2))
    cens=runif(n = N, min = 0, max = L)
    Obs=pmin(time,cens)
    Eps=(time<=cens)*Eps.fail
    New.Strata=ifelse(Eps==1, 1, ifelse(Z3==0, 2, 3))
    if(isFALSE(stratify)==T) Str_list=NULL else Str_list=stratified_sampling(str.var = New.Strata, size = c(length(which(New.Strata==1)), 75, 225))
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
      coverage[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(qnorm(Tau/P0), 0.5, -0.5+qnorm(Tau/P1)-qnorm(Tau/P0))
      tryCatch({
        NR=NewtonRhapson(COVAR = cbind(Z1,Z2), Obs = Obs, Eps = Eps, Tau = Tau, In.vec = in.vec, Str_list = Str_list)
        if(norm_vec(NR[[1]])>10){
          m.sol[[i]]=NULL
          boot.cov[[i]]=NULL
          coverage[[i]]=NULL
        } else {
          m.sol[[i]]=NR[[1]]
          boot.cov[[i]]=NR[[2]]
          coverage[[i]]=as.numeric(m.sol[[i]]-qnorm(p = 0.975)*sqrt(diag(boot.cov[[i]]))<in.vec & m.sol[[i]]+qnorm(p = 0.975)*sqrt(diag(boot.cov[[i]]))>in.vec)
          
        }
      }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    i=i+1
    if(length(Filter(Negate(is.null), m.sol))==m & length(Filter(Negate(is.null), boot.cov))==m) break
  }
  m.sol=Filter(Negate(is.null), m.sol)
  boot.cov=Filter(Negate(is.null), boot.cov)
  coverage=Filter(Negate(is.null), coverage)
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(m.sol)-1)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(boot.cov)*Reduce('+', boot.cov)
  return(list((1/length(m.sol)*Reduce('+', m.sol)-in.vec), sample.cov.mat, boot.mat, 1/length(coverage)*Reduce('+', coverage)))
}

NewtonRhapson=function(COVAR, Obs, Eps, Tau, In.vec=rep(0,ncol(COVAR)+1), Str_list=NULL){
  N=nrow(COVAR) ; P=ncol(COVAR)+1
  COVAR_int=cbind(rep(1, N), as.matrix(COVAR))

  repeat.beta=list(In.vec) ; repeat.cov=list(diag(rep(1/N, P))) ; j=2
  repeat{
    repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(obs = Obs, status = Eps, covariate = COVAR_int, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]), smooth.est.eq(obs = Obs, status = Eps, covariate = COVAR_int, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], tau = Tau, other.weight = Str_list[[2]])))
    V.cov=smooth.gamma.stratified(obs = Obs, status = Eps, covariate = COVAR_int, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], tau = Tau, strata_list = Str_list)
    repeat.cov[[j]]=1/N*solve(A(obs = Obs, status = Eps, covariate = COVAR_int, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]), V.cov)%*%solve(A(obs = Obs, status = Eps, covariate = COVAR_int, beta = repeat.beta[[j]], sigma = repeat.cov[[j-1]], other.weight = Str_list[[2]]))
    
    if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-3)|length(repeat.beta)==50) break else j=j+1
  }
  return(list(repeat.beta[[length(repeat.beta)]], repeat.cov[[length(repeat.cov)]]))
}

##################NWTSCO data##################
data(nwtsco, package = 'addhazard')
Time=c() ; Status=c()
for(i in 1:nrow(nwtsco)){
  if(nwtsco$relaps[i]==1 & nwtsco$dead[i]==0){Time[i]=nwtsco$trel[i] ; Status[i]=1}
  else if(nwtsco$relaps[i]==1 & nwtsco$dead[i]==1 & nwtsco$trel[i]<nwtsco$tsur[i]){Time[i]=nwtsco$trel[i] ; Status[i]=1}
  else if(nwtsco$relaps[i]==1 & nwtsco$dead[i]==1 & nwtsco$trel[i]==nwtsco$tsur[i]){Time[i]=nwtsco$tsur[i] ; Status[i]=1} 
  else {Time[i]=min(nwtsco$trel[i], nwtsco$tsur[i]) ; Status[i]=0}
}

Covariate=cbind(nwtsco$histol, nwtsco$age, ifelse(nwtsco$stage==2, 1, 0), ifelse(nwtsco$stage==3, 1, 0), ifelse(nwtsco$stage==4, 1, 0), ifelse(nwtsco$study==3, 0, 1))

Full_Cohort_005=NewtonRhapson(COVAR = Covariate, Obs = Time, Eps = Status, Tau = .05)
Full_Cohort_0075=NewtonRhapson(COVAR = Covariate, Obs = Time, Eps = Status, Tau = .075)
Full_Cohort_01=NewtonRhapson(COVAR = Covariate, Obs = Time, Eps = Status, Tau = .1)
Full_Cohort_0125=NewtonRhapson(COVAR = Covariate, Obs = Time, Eps = Status, Tau = .125)

Stage=ifelse(nwtsco$stage==1|nwtsco$stage==2, 0, 1)
Strata=c()
for(i in 1:nrow(nwtsco)){
  if(Status[i]==1) Strata[i]=0
  else if(Status[i]!=1 & Stage[i]==0 & nwtsco$age[i]<1 & nwtsco$histol[i]==0) Strata[i]=1
  else if(Status[i]!=1 & Stage[i]==0 & nwtsco$age[i]>=1 & nwtsco$histol[i]==0) Strata[i]=2
  else if(Status[i]!=1 & Stage[i]==1 & nwtsco$age[i]<1 & nwtsco$histol[i]==0) Strata[i]=3
  else if(Status[i]!=1 & Stage[i]==1 & nwtsco$age[i]>=1 & nwtsco$histol[i]==0) Strata[i]=4
  else if(Status[i]!=1 & Stage[i]==0 & nwtsco$age[i]<1 & nwtsco$histol[i]==1) Strata[i]=5
  else if(Status[i]!=1 & Stage[i]==0 & nwtsco$age[i]>=1 & nwtsco$histol[i]==1) Strata[i]=6
  else if(Status[i]!=1 & Stage[i]==1 & nwtsco$age[i]<1 & nwtsco$histol[i]==1) Strata[i]=7
  else if(Status[i]!=1 & Stage[i]==1 & nwtsco$age[i]>=1 & nwtsco$histol[i]==1) Strata[i]=8
}

Case_Cohort=function(rep, covariate, time, status, tau, strata){
  beta_list=list()
  Full_cohort_Var=NewtonRhapson(COVAR = covariate, Obs = time, Eps = status, Tau = tau)[[2]]
  j=1
  repeat{
    Strata_Sampled=stratified_sampling(str.var = strata, size = c(669, 120, 160, 25, 120, 12, 146, 5, 82))
    tryCatch({
      beta_list[[j]]=NewtonRhapson(COVAR = covariate, Obs = time, Eps = status, Tau = tau, In.vec = rep(0,ncol(covariate)+1), Str_list = Strata_Sampled)[[1]]
    }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
    j=j+1
    if(length(Filter(Negate(is.null), beta_list))==rep) break
  }
  beta_bar=as.vector(1/length(beta_list)*Reduce('+', beta_list))
  covariance_matrix_sum=lapply(X = beta_list, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/(length(beta_list)-1)*Reduce('+', covariance_matrix_sum)
  return(list(1/length(beta_list)*Reduce('+', beta_list), sample.cov.mat+Full_cohort_Var))
}

Case_Cohort_005=Case_Cohort(rep = 300, covariate = Covariate, time = Time, status = Status, tau = 0.05, strata = Strata)
Case_cohort_0075=Case_Cohort(rep = 300, covariate = Covariate, time = Time, status = Status, tau = 0.075, strata = Strata)
Case_Cohort_01=Case_Cohort(rep = 300, covariate = Covariate, time = Time, status = Status, tau = 0.1, strata = Strata)

##################WIHS data##################
data(wihs, package = "randomForestSRC")
mid_value=2
wihs$time_new=ifelse(wihs$time<=mid_value, wihs$time, mid_value)
wihs$time_new=wihs$time_new+runif(n = nrow(wihs), min = 0, max = 1.99999e-03)
wihs$status_new=ifelse(wihs$time<=mid_value, wihs$status, 0)

Covariate_wihs=wihs[,3:6]
Time_wihs=wihs$time_new
Status_wihs=wihs$status_new

Full_Cohort_wihs_0.1=NewtonRhapson(COVAR = Covariate_wihs, Obs = Time_wihs, Eps = Status_wihs, Tau = .1)
Full_Cohort_wihs_0.2=NewtonRhapson(COVAR = Covariate_wihs, Obs = Time_wihs, Eps = Status_wihs, Tau = .2)

##################Hodgkin's Disease data##################
data(hd, package = 'randomForestSRC')
hd$Sex=ifelse(hd$sex=='F', 1, 0)
hd$Treat=ifelse(hd$trtgiven=='RT', 0, 1)
hd$Extra=ifelse(hd$extranod=='Y', 1, 0)
hd$Clinic=hd$clinstg-1

attach(hd)

Full_Cohort_hd_005=NewtonRhapson(COVAR = cbind(age,Sex,Treat,Extra,Clinic), Obs = time, Eps = status, Tau = .05)
Full_Cohort_hd_01=NewtonRhapson(COVAR = cbind(age,Sex,Treat,Extra,Clinic), Obs = time, Eps = status, Tau = .1)
Full_Cohort_hd_015=NewtonRhapson(COVAR = cbind(age,Sex,Treat,Extra,Clinic), Obs = time, Eps = status, Tau = .15)
Full_Cohort_hd_02=NewtonRhapson(COVAR = cbind(age,Sex,Treat,Extra,Clinic), Obs = time, Eps = status, Tau = .2)
Full_Cohort_hd_025=NewtonRhapson(COVAR = cbind(age,Sex,Treat,Extra,Clinic), Obs = time, Eps = status, Tau = .25)
Full_Cohort_hd_03=NewtonRhapson(COVAR = cbind(age,Sex,Treat,Extra,Clinic), Obs = time, Eps = status, Tau = .3)

detach(hd)