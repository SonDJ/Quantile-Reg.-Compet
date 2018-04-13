ISMB_simulation_function=function(m, B, N, Dist, L, Tau){
  cov.est=diag(rep(1/N, 3), 3) ; m.sol=list() ; boot.cov=list() ; i=1
  repeat{
    Z1=runif(n = N, min = -1, max = 1) ; Z2=rbinom(n = N, size = 1, prob = 0.5)
    Eps=cause_sampling_func(Z2)
    time=sampling_func(dist = Dist, status = Eps, z1 = Z1, z2 = Z2)
    cens=runif(n = N, min = 0, max = L)
    Obs=ifelse(time>cens, cens, time)
    Eps=ifelse(Obs==cens, 0, Eps)
    
    if(length(survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)$time)<N) next
    
    else if(length(survfit(formula = Surv(time = Obs, event = ifelse(Eps==0, 1, 0))~1)$time)==N){
      pfcmp=suppressMessages(crrQR(ftime = log(Obs), fstatus = Eps, X = model.matrix(~Z1+Z2)[,-1], tau.range = c(Tau, Tau)))
      m.sol[[i]]=nleqslv(x = as.vector(pfcmp$beta.seq), fn = smooth.est.eq, method = c("Newton"), tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), sigma = cov.est)$x
      
      boot.list=list()
      
      for(j in 1:B){
        Eta=rexp(n = N, rate = 1)
        boot.list[[j]]=boot.smooth.est.eq(tau = Tau, n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), beta = m.sol[[i]], sigma = cov.est, eta = Eta)
      }
      
      boot.bar=as.vector(1/length(boot.list)*Reduce('+', boot.list))
      boot_matrix_sum=lapply(X = boot.list, FUN = function(j){outer(j-boot.bar, j-boot.bar)})
      V.cov=1/B*Reduce('+', boot_matrix_sum)
      A_mat=A(n = N, obs = Obs, status = Eps, covariate = cbind(rep(1, N), Z1, Z2), beta = m.sol[[i]], sigma = cov.est)
      boot.cov[[i]]=solve(A_mat)%*%V.cov%*%solve(A_mat)
    }
    i=i+1
    if(length(m.sol)==m & length(boot.cov)==m) break
  }
  beta_bar=as.vector(1/length(m.sol)*Reduce('+', m.sol))
  covariance_matrix_sum=lapply(X = m.sol, FUN = function(j){outer(j-beta_bar, j-beta_bar)})
  sample.cov.mat=1/length(m.sol)*Reduce('+', covariance_matrix_sum)
  boot.mat=1/length(m.sol)*Reduce('+', boot.cov)
  return(list(1/length(m.sol)*Reduce('+', m.sol), sample.cov.mat, boot.mat))
}