library(survival)
na.omit.list=function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
inverse.check=function(m) class(try(solve(m),silent=T))=="matrix"
quadform=function(A,x) {return(colSums(x * (A %*% x)))}
norm_vec=function(x) return(sqrt(sum(x^2)))
quadform.vec=function(Z, sig) {return(diag(Z%*%sig%*%t(Z)))}
censoring_rate=function(N, Dist, Cohort, L, numrep=1e+03){
  num.cens=c()
  for(i in 1:numrep){
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2, p0 = ifelse(Cohort==T, 0.9, 0.7), p1 = ifelse(Cohort==T, 0.7, 0.8))
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    if(Cohort==T){
      cens=mixture_censoring_sample(num = N, l = L)
    } else{
      cens=runif(n = N, min = 0, max = L)
    }
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

mixture_censoring_sample=function(num, l){
  unif=runif(n = num, min = 0, max = 1)
  censoring=c()
  for(i in 1:num){
    if(unif[i]>=0.8) censoring[i]=l
    else censoring[i]=l*unif[i]
  }
  return(censoring)
}

A=function(n, obs, status, covariate, beta, sigma){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) ; obs=obs[obs.per]
  status=status[obs.per] ; covariate=covariate[obs.per,]
  
  no.na=setdiff(x = 1:n, which(is.na(covariate[,1])))
  obs=obs[no.na]
  status=status[no.na]
  km.cens=km.cens[no.na]
  covariate=covariate[no.na,]
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  
  dnorm.vec=dnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)/sqrt(quadvec)
  weight=ifelse(status==1, 1/km.cens, 0)
  return(1/n*t(covariate)%*%diag(weight*dnorm.vec)%*%covariate)
}

smooth.gamma=function(beta, tau, n, obs, status, covariate, sigma, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) ; obs=obs[obs.per]
  status=status[obs.per] ; covariate=covariate[obs.per,]
  
  no.na=setdiff(x = 1:n, which(is.na(covariate[,1])))
  obs=obs[no.na]
  status=status[no.na]
  km.cens=km.cens[no.na]
  covariate=covariate[no.na,]
  if(is.null(CCH)==T) pn=1 else pn=CCH
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  cc.weight=ifelse(status==1, 1, 1/pn)
  
  mat=sweep(covariate, MARGIN = 1, cc.weight*(weight.vec*pnorm.vec-rep(tau, length(status))), '*')
  
  cen.vec=ifelse(status==0, 1, 0)
  l=matrix(rep(0, length(status)*ncol(covariate)), nrow = length(status))
  for(i in 1:length(status)){
    l[i,]=colSums(sweep(covariate, MARGIN = 1, ifelse(obs>=obs[i], 1, 0)*pnorm.vec*weight.vec, FUN = '*'))/sum(ifelse(obs>=obs[i], 1, 0))
  }

  return(1/(n^2)*t(mat)%*%mat-1/(n^2)*t(l)%*%diag(cen.vec)%*%l)
}

smooth.gamma.new=function(beta, tau, n, obs, status, covariate, sigma, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) ; obs=obs[obs.per]
  status=status[obs.per] ; covariate=covariate[obs.per,]
  
  no.na=setdiff(x = 1:n, which(is.na(covariate[,1])))
  obs=obs[no.na]
  status=status[no.na]
  km.cens=km.cens[no.na]
  covariate=covariate[no.na,]
  if(is.null(CCH)==T) pn=1 else pn=CCH
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  cc.weight=ifelse(status==1, 1, 1/pn)
  
  mat=sweep(covariate, MARGIN = 1, (weight.vec*pnorm.vec-rep(tau, length(status))), '*')
  
  cen.vec=ifelse(status==0, 1, 0)
  l=matrix(rep(0, length(status)*ncol(covariate)), nrow = length(status))
  for(i in 1:length(status)){
    l[i,]=colSums(sweep(covariate, MARGIN = 1, ifelse(obs>=obs[i], 1, 0)*pnorm.vec*weight.vec, FUN = '*'))/sum(ifelse(obs>=obs[i], 1, 0))
  }
  
  return(1/(n^2)*t(mat-l)%*%diag(cc.weight)%*%(mat-l)+(1/n^2)*((1-pn)/pn)*t(mat)%*%diag(ifelse(status==1, 0, 1/pn))%*%mat)
}

smooth.est.eq=function(beta, tau, n, obs, status, covariate, sigma, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) ; obs=obs[obs.per]
  status=status[obs.per] ; covariate=covariate[obs.per,]
  
  no.na=setdiff(x = 1:n, which(is.na(covariate[,1])))
  obs=obs[no.na]
  status=status[no.na]
  km.cens=km.cens[no.na]
  covariate=covariate[no.na,]
  if(is.null(CCH)==T) pn=1 else pn=CCH
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  cc.weight=ifelse(status==1, 1, 1/pn)
  
  return(1/n*colSums(sweep(covariate, MARGIN = 1, cc.weight*(weight.vec*pnorm.vec-rep(tau, length(status))), '*')))
}

#MB method
boot.smooth.est.eq=function(tau, n, obs, status, covariate, beta, sigma, eta, CCH=NULL){
  SF=survfit(formula = Surv(time = obs, event = ifelse(status==0, 1, 0))~1)
  km.cens=SF$surv
  obs.per=order(obs) ; obs=obs[obs.per]
  status=status[obs.per] ; covariate=covariate[obs.per,]
  
  no.na=setdiff(x = 1:n, which(is.na(covariate[,1])))
  obs=obs[no.na]
  status=status[no.na]
  km.cens=km.cens[no.na]
  covariate=covariate[no.na,]
  eta=eta[no.na]
  if(is.null(CCH)==T) pn=1 else pn=CCH
  
  quadvec=quadform.vec(Z = covariate, sig = sigma)
  weight.vec=ifelse(status==1, 1/km.cens, 0)
  pnorm.vec=pnorm(-(log(obs)-as.vector(covariate%*%beta))/sqrt(quadvec), 0, 1)
  cc.weight=ifelse(status==1, 1, 1/pn)
  
  return(1/n*colSums(sweep(covariate, MARGIN = 1, eta*cc.weight*(weight.vec*pnorm.vec-rep(tau, length(status))), '*')))
}

#Case-cohort design
CCHD=function(rate, data, status){
  location=which(status==1)
  sub.coh.ind=sample(x = 1:nrow(data), size = floor(nrow(data)*rate), replace = F)
  case.coh.ind=union(location, sub.coh.ind)
  data[setdiff(x = 1:nrow(data), y = case.coh.ind), ]=rep(NA, ncol(data))
  return(list(data, rate, sub.coh.ind))
}

# Iterative Method
Iter_simulation_gamma=function(m, N, Dist, L, Tau, Cohort=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    P0=0.7 ; P1=0.8
    Eps.fail=cause_sampling_func(Z2, p0 = P0, p1 = P1)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    CCH_desig=CCHD(rate = 0.25, data = cbind(rep(1, N), Z1, Z2), status = Eps)
    if(Cohort==T) COVAR=CCH_desig[[1]] else COVAR=cbind(rep(1, N), Z1, Z2)
    if(Cohort==T) CCH_status=0.25 else CCH_status=NULL
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      in.vec=c(-1+qnorm(Tau/P0), 1, 1+qnorm(Tau/P1)-qnorm(Tau/P0))
      repeat.beta=list(in.vec) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]))==F|any(is.nan(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]])))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]), smooth.est.eq(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], CCH = CCH_status, tau = Tau)))
          
          V.cov=smooth.gamma.new(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], CCH = CCH_status, tau = Tau)
          
          repeat.cov[[j]]=solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]), V.cov)%*%solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]))
          if(norm_vec(repeat.beta[[j]]-repeat.beta[[j-1]])<10^(-4)|length(repeat.beta)==50) break else j=j+1
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

Iter_simulation=function(m, B, N, Dist, L, Tau, Cohort=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    CCH_desig=CCHD(rate = 0.1, data = cbind(rep(1, N), Z1, Z2), status = Eps)
    if(Cohort==T) COVAR=CCH_desig[[1]] else COVAR=cbind(rep(1, N), Z1, Z2)
    if(Cohort==T) CCH_status=0.1 else CCH_status=NULL
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      repeat.beta=list(c(-1+qnorm(Tau/0.7), 1, 1+qnorm(Tau/0.8)-qnorm(Tau/0.7))) ; repeat.cov=list(naive.cov) ; j=2
      repeat{
        if(inverse.check(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]))==F|any(is.nan(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]])))==T) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]), smooth.est.eq(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], CCH = CCH_status, tau = Tau)))
          
          boot.list=list()
          for(k in 1:B){
            Eta=rexp(n = N, rate = 1)
            boot.list[[k]]=boot.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], eta = Eta)
          }
          
          boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
          boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
          
          V.cov=1/(B-1)*Reduce('+', boot_matrix_sum)          
          repeat.cov[[j]]=solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]), V.cov)%*%solve(A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]]))
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

Iter_simulation_A=function(m, B, N, Dist, L, Tau, Cohort=F){
  naive.cov=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>=time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    CCH_desig=CCHD(rate = 0.1, data = cbind(rep(1, N), Z1, Z2), status = Eps)
    if(Cohort==T) COVAR=CCH_desig[[1]] else COVAR=cbind(rep(1, N), Z1, Z2)
    if(Cohort==T) CCH_status=CCH_desig[[2]] else CCH_status=NULL
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    Obs.per=match(x = KME$time, table = Obs)
    Weight=ifelse(Eps[Obs.per]==1, 1/KME$surv, 0)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      repeat.beta=list(c(-1+qnorm(Tau/0.7), 1, 1+qnorm(Tau/0.8)-qnorm(Tau/0.7))) ; repeat.cov=list(naive.cov) ; j=2 ; A.row=list()
      repeat{
        Z=mvrnorm(n = B, mu = rep(0, ncol(COVAR)), Sigma = diag(rep(1, ncol(COVAR))))
        A=matrix(data = rep(0, ncol(COVAR)^2), ncol(COVAR))
        for(l in 1:ncol(COVAR)){
          y=c()
          for(t in 1:B){
            y[t]=(N^(1/2))*smooth.est.eq(beta = (repeat.beta[[j-1]]+(N^(-1/2))*Z[t,]), tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, sigma = repeat.cov[[j-1]], CCH = CCH_status)[l]
          }
          df = cbind(y, Z) ; df=as.data.frame(df)
          colnames(df)=c('y', paste0('x',1:ncol(COVAR)))
          A[l,]=as.vector(lm(formula = y~x1+x2+x3+0, data = df)$coefficients)
        }
        
        if(inverse.check(A)==F) {repeat.beta[[j]]=NA ; break}
        
        else{
          repeat.beta[[j]]=repeat.beta[[j-1]]-as.vector(solve(A)%*%smooth.est.eq(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], CCH = CCH_status, tau = Tau))
          boot.list=list()
          for(k in 1:B){
            Eta=rexp(n = N, rate = 1)
            boot.list[[k]]=boot.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, beta = repeat.beta[[j-1]], sigma = repeat.cov[[j-1]], eta = Eta)
          }
          
          boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
          boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
          
          V.cov=1/(B-1)*Reduce('+', boot_matrix_sum)
          repeat.cov[[j]]=solve(A)%*%V.cov%*%solve(A)
          if(max(abs((repeat.beta[[j]]-repeat.beta[[j-1]])))<10^(-2)|length(repeat.beta)==50) break else j=j+1
        }
      }
      
      if(any(is.na(repeat.beta[[length(repeat.beta)]]))==T) {m.sol[[i]]=NULL ; boot.cov[[i]]=NULL}
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

#IS-MB method
ISMB_simulation_function=function(m, B, N, Dist, L, Tau, Cohort=F){
  cov.est=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps.fail=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps.fail, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(cens>=time, time, cens)
    Eps=ifelse(Obs==cens, 0, Eps.fail)
    CCH_desig=CCHD(rate = 0.2, data = cbind(rep(1, N), Z1, Z2), status = Eps)
    COVAR=if(Cohort==T) CCH_desig[[1]] else cbind(rep(1, N), Z1, Z2)
    if(Cohort==T) CCH_status=CCH_desig[[2]] else CCH_status=NULL
    KME=survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)
    Obs.per=match(x = KME$time, table = Obs)
    Weight=ifelse(Eps[Obs.per]==1, 1/KME$surv, 0)
    
    if(length(KME$time)<N){
      m.sol[[i]]=NULL
      boot.cov[[i]]=NULL
    }
    
    else if(length(KME$time)==N){
      
      m.sol[[i]]=nleqslv(x = c(-1.431, 1, 0.756), fn = smooth.est.eq, tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, sigma = cov.est, CCH = CCH_status)$x
      
      boot.list=list() ; Eta=list()
      
      for(j in 1:B){
        Eta[[j]]=rexp(n = N, rate = 1)
        boot.list[[j]]=boot.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = COVAR, beta = m.sol[[i]], sigma = cov.est, eta = Eta[[j]], CCH = CCH_status)
      }
      
      boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
      boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
      V.cov=1/(B-1)*Reduce('+', boot_matrix_sum)
      A_mat=A(n = N, obs = Obs, status = Eps, covariate = COVAR, beta = m.sol[[i]], sigma = cov.est, CCH = CCH_status)
      if(rcond(A_mat)>=1e-15) solve.A = solve(A_mat) else solve.A=ginv(A_mat)
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
