
rm(list = ls())

# setwd("~/TE/s1")
# n.itt = 10000
# n.B = 1000
# n.cluster = 120
# scen.all.ind = 1

setwd("~/Dropbox/Research/AbbVie/Treatment_effect/R code/SIM_R1_code/2_Sim_cond/")
# setwd("C:/Users/ZHANTX/Documents/TZ/code/treatment effect/R code/Share/")
n.itt = 10^4
n.B = 80
n.cluster = 8

library(doParallel)

n = 40
time.start = Sys.time()
n.scen = 4

final.out.mat = matrix(NA, nrow = n.scen, ncol = 25)

for (scen.all.ind in 1:4){
for (scen.itt in 1:n.scen){

if (scen.all.ind==1){ 
  gamma.type = "NA"
if (scen.itt==1) {mu.true.vec = c(1, 1, 1); sd.true.vec = c(5, 5, 5); gamma.prop = 0 }
if (scen.itt==2) {mu.true.vec = c(1, 1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0 }
if (scen.itt==3) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0 }
if (scen.itt==4) {mu.true.vec = c(1, 1.2, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0 }
}
  
if (scen.all.ind==2){ 
  gamma.type = "NA"
if (scen.itt==1) {mu.true.vec = c(1, 1, 1); sd.true.vec = c(3, 4, 5); gamma.prop = 0 }
if (scen.itt==2) {mu.true.vec = c(1, 1, 1.2); sd.true.vec = c(3, 4, 5) ; gamma.prop = 0}
if (scen.itt==3) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(3, 4, 5); gamma.prop = 0 }
if (scen.itt==4) {mu.true.vec = c(1, 1.2, 1.2); sd.true.vec = c(3, 4, 5); gamma.prop = 0 }
}
    
if (scen.all.ind==3){
  gamma.type = "gamma"
if (scen.itt==1) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.1}
if (scen.itt==2) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.2}
if (scen.itt==3) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.3} 
if (scen.itt==4) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.5}
}
  
if (scen.all.ind==4){
  gamma.type = "unif"
  if (scen.itt==1) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.1}
  if (scen.itt==2) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.2}
  if (scen.itt==3) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.3} 
  if (scen.itt==4) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5); gamma.prop = 0.5}
}

gamma.scale.vec = sd.true.vec^2/mu.true.vec
gamma.shape.vec = mu.true.vec/gamma.scale.vec

unif.b.vec = (sqrt(12)*sd.true.vec + 2*mu.true.vec)/2
unif.a.vec = (-sqrt(12)*sd.true.vec + 2*mu.true.vec)/2
  
# a = rgamma(10000, shape = gamma.shape.vec[1], scale = gamma.scale.vec[1])
# b = rnorm(100000, mean = mu.true.vec[1], sd = sd.true.vec[1])
# 
# mean(a)
# sd(a)

############################## Parametric bootstrap
cl = makeCluster(n.cluster)
registerDoParallel(cl)
out.vec = foreach(itt=1:n.itt) %dopar% {  
  
  source("sim_cond_func.R")
  set.seed(itt+n.itt*scen.itt)
  
  if (scen.all.ind<=3){
  mu1.data = (1-gamma.prop)*rnorm(n, mu.true.vec[1], sd = sd.true.vec[1])+
    gamma.prop*rgamma(n, shape = gamma.shape.vec[1], scale = gamma.scale.vec[1])
  mu2.data = (1-gamma.prop)*rnorm(n, mu.true.vec[2], sd = sd.true.vec[2])+
    gamma.prop*rgamma(n, shape = gamma.shape.vec[2], scale = gamma.scale.vec[2])
  mu3.data = (1-gamma.prop)*rnorm(n, mu.true.vec[3], sd = sd.true.vec[3])+
    gamma.prop*rgamma(n, shape = gamma.shape.vec[3], scale = gamma.scale.vec[3])
  } else{
  mu1.data = (1-gamma.prop)*rnorm(n, mu.true.vec[1], sd = sd.true.vec[1])+
    gamma.prop*runif(n, min = unif.a.vec[1], max = unif.b.vec[1])
  mu2.data = (1-gamma.prop)*rnorm(n, mu.true.vec[2], sd = sd.true.vec[2])+
    gamma.prop*runif(n, min = unif.a.vec[2], max = unif.b.vec[2])
  mu3.data = (1-gamma.prop)*rnorm(n, mu.true.vec[3], sd = sd.true.vec[3])+
    gamma.prop*runif(n, min = unif.a.vec[3], max = unif.b.vec[3])
  }
  
  mu1.sample = mean(mu1.data)
  mu2.sample = mean(mu2.data)
  mu3.sample = mean(mu3.data)
  
  if ((mu3.sample<mu2.sample)|(mu3.sample<mu1.sample)){
    PB.est.1.out = PB.est.2.out = PB.S.est.2.out = 
      NB.est.1.out = NB.est.2.out = NB.S.est.2.out = 
      JK.est.out = est.L.out  = NA
  } else {
    PB.est.1.out = PB.est.1.func(mu1.data, mu2.data, mu3.data, n, n.B)
    PB.est.2.out = PB.est.2.func(mu1.data, mu2.data, mu3.data, n, n.B)
    PB.S.est.2.out = PB.S.est.2.func(mu1.data, mu2.data, mu3.data, n, n.B)
    
    NB.est.1.out = NB.est.1.func(mu1.data, mu2.data, mu3.data, n, n.B)
    NB.est.2.out = NB.est.2.func(mu1.data, mu2.data, mu3.data, n, n.B)
    NB.S.est.2.out = NB.S.est.2.func(mu1.data, mu2.data, mu3.data, n, n.B)
    
    JK.est.out = JK.est.func(mu1.data, mu2.data, mu3.data, n)
    
    est.L.out = est.L.func(mu1.data, mu2.data, mu3.data, n)
  }
  
  return(c(mu1.sample, 
           mu2.sample,
           mu3.sample, 
           
           max(mu1.sample, mu2.sample, mu3.sample),
           
           PB.est.1.out,
           PB.est.2.out, 
           PB.S.est.2.out, 
           NB.est.1.out,
           NB.est.2.out, 
           NB.S.est.2.out,
           
           JK.est.out,
           est.L.out
  ))
}
stopCluster(cl)

out.mat.temp = matrix(unlist(out.vec),nrow = n.itt, ncol=12, byrow = TRUE)
out.mat = out.mat.temp[!is.na(out.mat.temp[,6]), ]

out.final.vec = c(
                     ## bias
                     apply(out.mat, 2, mean)-c(mu.true.vec, 
                                   rep(max(mu.true.vec), 9)),
                     ## MSE
                     mean((out.mat[,4]-max(mu.true.vec))^2),
                     mean((out.mat[,5]-max(mu.true.vec))^2),
                     mean((out.mat[,6]-max(mu.true.vec))^2),
                     mean((out.mat[,7]-max(mu.true.vec))^2),
                     mean((out.mat[,8]-max(mu.true.vec))^2),
                     mean((out.mat[,9]-max(mu.true.vec))^2),
                     mean((out.mat[,10]-max(mu.true.vec))^2),
                     mean((out.mat[,11]-max(mu.true.vec))^2),
                     mean((out.mat[,12]-max(mu.true.vec))^2))
                     
                     
print(out.final.vec)

final.out.mat[scen.itt, ] = c(paste(as.character(mu.true.vec), collapse = ", "),
                              paste(as.character(sd.true.vec), collapse = ", "),
                              gamma.prop, 
                              gamma.type,
                              out.final.vec
                              )
colnames(final.out.mat) = c("mu", "sd", "gamma_prop", "gamma_type",
                            "bias_1", "bias_2", "bias_3", 
                            "bias_naive", 
                            "bias_PB_S", "bias_PB_D", "bias_PB_D_S",
                            "bias_NB_S", "bias_NB_D", "bias_NB_D_S",
                            "bias_JK", "bias_L", 
                            "MSE_naive", 
                            "MSE_PB_S", "MSE_PB_D", "MSE_PB_D_S",
                            "MSE_NB_S", "MSE_NB_D", "MSE_NB_D_S",
                            "MSE_JK", "MSE_L"
)
}

print(Sys.time()-time.start)
write.csv(final.out.mat, paste0(scen.all.ind, "_results_cond_", n.itt, ".csv"))
}





