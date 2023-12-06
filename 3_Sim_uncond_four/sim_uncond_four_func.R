
## parametric bootstrap single
PB.est.1.func = function(mu1.data.1.in, mu2.data.1.in, mu3.data.1.in, 
                         mu4.data.1.in, 
                      n.in, n.B.in){
  mu12.boot = sapply(1:n.B.in, function(x1){
    mu1.boot.1 = rnorm(n.in, mean(mu1.data.1.in), sd = sd(mu1.data.1.in))
    mu2.boot.1 = rnorm(n.in, mean(mu2.data.1.in), sd = sd(mu2.data.1.in)) 
    mu3.boot.1 = rnorm(n.in, mean(mu3.data.1.in), sd = sd(mu3.data.1.in)) 
    mu4.boot.1 = rnorm(n.in, mean(mu4.data.1.in), sd = sd(mu4.data.1.in)) 
    
    return(max(mean(mu1.boot.1), mean(mu2.boot.1), mean(mu3.boot.1),
               mean(mu4.boot.1)))
  })
  return(2*max(mean(mu1.data.1.in), mean(mu2.data.1.in), mean(mu3.data.1.in),
               mean(mu4.data.1.in))-
           mean(mu12.boot))
}

## parametric bootstrap double
PB.est.2.func = function(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                         mu4.data.2.in,
                      n.in, n.B.in){
  mu12.boot.2 = sapply(1:n.B.in, function(x2){
    mu1.boot.2 = rnorm(n.in, mean(mu1.data.2.in), sd = sd(mu1.data.2.in))
    mu2.boot.2 = rnorm(n.in, mean(mu2.data.2.in), sd = sd(mu2.data.2.in)) 
    mu3.boot.2 = rnorm(n.in, mean(mu3.data.2.in), sd = sd(mu3.data.2.in)) 
    mu4.boot.2 = rnorm(n.in, mean(mu4.data.2.in), sd = sd(mu4.data.2.in)) 
    
    return(PB.est.1.func(mu1.boot.2, mu2.boot.2, mu3.boot.2, 
                         mu4.boot.2, n.in, n.B.in))
  })
  return(2*PB.est.1.func(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                         mu4.data.2.in, 
                      n.in, n.B.in)-
           mean(mu12.boot.2))
}

## parametric bootstrap double with shrinkage
PB.S.est.2.func = function(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                         mu4.data.2.in,
                         n.in, n.B.in){
  mu12.boot.2 = sapply(1:n.B.in, function(x2){
    mu1.boot.2 = rnorm(n.in, mean(mu1.data.2.in), sd = sd(mu1.data.2.in))
    mu2.boot.2 = rnorm(n.in, mean(mu2.data.2.in), sd = sd(mu2.data.2.in)) 
    mu3.boot.2 = rnorm(n.in, mean(mu3.data.2.in), sd = sd(mu3.data.2.in)) 
    mu4.boot.2 = rnorm(n.in, mean(mu4.data.2.in), sd = sd(mu4.data.2.in)) 
    
    return(PB.est.1.func(mu1.boot.2, mu2.boot.2, mu3.boot.2, 
                         mu4.boot.2, n.in, n.B.in))
  })
  
  XS.func = 2*PB.est.1.func(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                            mu4.data.2.in, 
                            n.in, n.B.in)-
    mean(mu12.boot.2)
  Xbar.func = mean(c(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                     mu4.data.2.in))
  sigma2.func = mean(var(mu1.data.2.in), var(mu2.data.2.in),
                     var(mu3.data.2.in), var(mu4.data.2.in))
  C.func = max(0, 1-3*sigma2.func/(n.in*((mean(mu1.data.2.in)-Xbar.func)^2+
                                           (mean(mu2.data.2.in)-Xbar.func)^2+
                                           (mean(mu3.data.2.in)-Xbar.func)^2+
                                           (mean(mu4.data.2.in)-Xbar.func)^2
  )))
  
  return(C.func*XS.func + (1-C.func)*Xbar.func)

}

## non-parametric bootstrap single
NB.est.1.func = function(mu1.data.1.in, mu2.data.1.in, mu3.data.1.in,
                         mu4.data.1.in, 
                         n.in, n.B.in){
  mu12.boot = sapply(1:n.B.in, function(x1){
    mu1.boot.1 = sample(mu1.data.1.in, n.in, replace = TRUE)
    mu2.boot.1 = sample(mu2.data.1.in, n.in, replace = TRUE)
    mu3.boot.1 = sample(mu3.data.1.in, n.in, replace = TRUE)
    mu4.boot.1 = sample(mu4.data.1.in, n.in, replace = TRUE)
    
    return(max(mean(mu1.boot.1), mean(mu2.boot.1), mean(mu3.boot.1),
               mean(mu4.boot.1)))
  })
  return(2*max(mean(mu1.data.1.in), mean(mu2.data.1.in), mean(mu3.data.1.in),
               mean(mu4.data.1.in))-
           mean(mu12.boot))
}

## non-parametric bootstrap double
NB.est.2.func = function(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in, 
                         mu4.data.2.in, n.in, n.B.in){
  mu12.boot.2 = sapply(1:n.B.in, function(x2){
    mu1.boot.2 = sample(mu1.data.2.in, n.in, replace = TRUE)
    mu2.boot.2 = sample(mu2.data.2.in, n.in, replace = TRUE)
    mu3.boot.2 = sample(mu3.data.2.in, n.in, replace = TRUE)
    mu4.boot.2 = sample(mu4.data.2.in, n.in, replace = TRUE)
    
    return(NB.est.1.func(mu1.boot.2, mu2.boot.2, mu3.boot.2, 
                         mu4.boot.2, n.in, n.B.in))
  })
  return(2*NB.est.1.func(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                         mu4.data.2.in, 
                         n.in, n.B.in)-
           mean(mu12.boot.2))
}

## non-parametric bootstrap double with Shrinkage
NB.S.est.2.func = function(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in, 
                         mu4.data.2.in, n.in, n.B.in){
  mu12.boot.2 = sapply(1:n.B.in, function(x2){
    mu1.boot.2 = sample(mu1.data.2.in, n.in, replace = TRUE)
    mu2.boot.2 = sample(mu2.data.2.in, n.in, replace = TRUE)
    mu3.boot.2 = sample(mu3.data.2.in, n.in, replace = TRUE)
    mu4.boot.2 = sample(mu4.data.2.in, n.in, replace = TRUE)
    
    return(NB.est.1.func(mu1.boot.2, mu2.boot.2, mu3.boot.2, 
                         mu4.boot.2, n.in, n.B.in))
  })

  XS.func = 2*NB.est.1.func(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                            mu4.data.2.in, 
                            n.in, n.B.in)-
    mean(mu12.boot.2)
  Xbar.func = mean(c(mu1.data.2.in, mu2.data.2.in, mu3.data.2.in,
                     mu4.data.2.in))
  sigma2.func = mean(var(mu1.data.2.in), var(mu2.data.2.in),
                     var(mu3.data.2.in), var(mu4.data.2.in))
  C.func = max(0, 1-3*sigma2.func/(n.in*((mean(mu1.data.2.in)-Xbar.func)^2+
                                           (mean(mu2.data.2.in)-Xbar.func)^2+
                                           (mean(mu3.data.2.in)-Xbar.func)^2+
                                           (mean(mu4.data.2.in)-Xbar.func)^2
  )))
  
  return(C.func*XS.func + (1-C.func)*Xbar.func)
  
}

## Jackknife
JK.est.func = function(mu1.data.1.in, mu2.data.1.in, mu3.data.1.in, 
                       mu4.data.1.in, n.in){

  mu12.JK = sapply(1:n.in, function(x1.in){
    mu1.data.2.in = mu1.data.1.in[-x1.in]
    mu2.data.2.in = mu2.data.1.in[-x1.in]
    mu3.data.2.in = mu3.data.1.in[-x1.in]
    mu4.data.2.in = mu4.data.1.in[-x1.in]
    
    mu1.sample.2 = mean(mu1.data.2.in)
    mu2.sample.2 = mean(mu2.data.2.in)
    mu3.sample.2 = mean(mu3.data.2.in)
    mu4.sample.2 = mean(mu4.data.2.in)
    
    return(max(mu1.sample.2, mu2.sample.2, mu3.sample.2, mu4.sample.2))
  })
  
  jack.1 = n.in*max(mean(mu1.data.1.in), mean(mu2.data.1.in), 
                 mean(mu3.data.1.in), mean(mu4.data.1.in)) - 
    (n.in-1)*mean(mu12.JK)
  
  return(jack.1)
}

## Lindley's shrinkage estimator
est.L.func = function(mu1.data.L.in, mu2.data.L.in, mu3.data.L.in, 
                      mu4.data.L.in, n.in){
  XS.func = max(mean(mu1.data.L.in), mean(mu2.data.L.in), mean(mu3.data.L.in),
                mean(mu4.data.L.in))
  Xbar.func = mean(c(mu1.data.L.in, mu2.data.L.in, mu3.data.L.in,
                     mu4.data.L.in))
  sigma2.func = mean(var(mu1.data.L.in), var(mu2.data.L.in),
                     var(mu3.data.L.in), var(mu4.data.L.in))
  C.func = max(0, 1-3*sigma2.func/(n.in*((mean(mu1.data.L.in)-Xbar.func)^2+
                               (mean(mu2.data.L.in)-Xbar.func)^2+
                               (mean(mu3.data.L.in)-Xbar.func)^2+
                                 (mean(mu4.data.L.in)-Xbar.func)^2   
  )))
  # C.func = 0
  return(C.func*XS.func + (1-C.func)*Xbar.func)
}











