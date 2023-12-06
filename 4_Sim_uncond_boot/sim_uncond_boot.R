
rm(list = ls())

# setwd("~/TE/R1/")
# n.itt = 10000
# n.cluster = 120
# scen.all.ind = 1

setwd("~/Dropbox/Research/AbbVie/Treatment_effect/R code/SIM_R1_code/4_Sim_uncond_boot/")
# setwd("C:/Users/ZHANTX/Documents/TZ/code/treatment effect/R code/Share/")
n.itt = 10^4
n.cluster = 8
scen.all.ind = 1

library(doParallel)

n = 40
time.start = Sys.time()
n.scen = 4

gamma.type = "NA"
gamma.prop = 0 

if (scen.all.ind==1) {mu.true.vec = c(1, 1, 1); sd.true.vec = c(5, 5, 5)}
if (scen.all.ind==2) {mu.true.vec = c(1, 1, 1.2); sd.true.vec = c(5, 5, 5)}
if (scen.all.ind==3) {mu.true.vec = c(1, 1.1, 1.2); sd.true.vec = c(5, 5, 5)}
if (scen.all.ind==4) {mu.true.vec = c(1, 1.2, 1.2); sd.true.vec = c(5, 5, 5)}

final.out.mat = matrix(NA, nrow = n.scen, ncol = 22)

for (scen.itt in 1:n.scen){
print(scen.itt)

# if (scen.itt==1){n.B.S = 81^3; n.B.D = 729; n.B.T = 81}
# if (scen.itt==2){n.B.S = 10^6; n.B.D = 10^3; n.B.T = 100}

if (scen.itt==1){n.B.S = 80; n.B.D = 80; n.B.T = 80}
if (scen.itt==2){n.B.S = 100; n.B.D = 100; n.B.T = 100}

if (scen.itt==3){n.B.S = 500; n.B.D = 500; n.B.T = 2}
if (scen.itt==4){n.B.S = 1000; n.B.D = 1000; n.B.T = 2}

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
  
  source("sim_uncond_func.R")
  set.seed(itt+n.itt*scen.itt)
 
  mu1.data = (1-gamma.prop)*rnorm(n, mu.true.vec[1], sd = sd.true.vec[1])+
    gamma.prop*rgamma(n, shape = gamma.shape.vec[1], scale = gamma.scale.vec[1])
  mu2.data = (1-gamma.prop)*rnorm(n, mu.true.vec[2], sd = sd.true.vec[2])+
    gamma.prop*rgamma(n, shape = gamma.shape.vec[2], scale = gamma.scale.vec[2])
  mu3.data = (1-gamma.prop)*rnorm(n, mu.true.vec[3], sd = sd.true.vec[3])+
    gamma.prop*rgamma(n, shape = gamma.shape.vec[3], scale = gamma.scale.vec[3])
  
  mu1.sample = mean(mu1.data)
  mu2.sample = mean(mu2.data)
  mu3.sample = mean(mu3.data)
  
  PB.est.1.out = PB.est.1.func(mu1.data, mu2.data, mu3.data, n, n.B.S)
  PB.est.2.out = PB.est.2.func(mu1.data, mu2.data, mu3.data, n, n.B.D)
  PB.est.3.out = PB.est.3.func(mu1.data, mu2.data, mu3.data, n, n.B.T)
  
  NB.est.1.out = NB.est.1.func(mu1.data, mu2.data, mu3.data, n, n.B.S)
  NB.est.2.out = NB.est.2.func(mu1.data, mu2.data, mu3.data, n, n.B.D)
  NB.est.3.out = NB.est.3.func(mu1.data, mu2.data, mu3.data, n, n.B.T)

  return(c(mu1.sample, 
           mu2.sample,
           mu3.sample, 
           
           PB.est.1.out,
           PB.est.2.out, 
           PB.est.3.out, 
           NB.est.1.out,
           NB.est.2.out,
           NB.est.3.out
  ))
}
stopCluster(cl)

out.mat = matrix(unlist(out.vec),nrow = n.itt, ncol=9, byrow = TRUE)

out.final.vec = c(
                     ## bias
                     apply(out.mat, 2, mean)-c(mu.true.vec, 
                                   rep(max(mu.true.vec), 6)),
                     ## MSE
                     mean((out.mat[,4]-max(mu.true.vec))^2),
                     mean((out.mat[,5]-max(mu.true.vec))^2),
                     mean((out.mat[,6]-max(mu.true.vec))^2),
                     mean((out.mat[,7]-max(mu.true.vec))^2),
                     mean((out.mat[,8]-max(mu.true.vec))^2),
                     mean((out.mat[,9]-max(mu.true.vec))^2))
                     
                     
print(out.final.vec)

final.out.mat[scen.itt, ] = c(paste(as.character(mu.true.vec), collapse = ", "),
                              paste(as.character(sd.true.vec), collapse = ", "),
                              gamma.prop, 
                              gamma.type,
                              n.B.S,
                              n.B.D,
                              n.B.T,
                              out.final.vec
                              )
colnames(final.out.mat) = c("mu", "sd", "gamma_prop", "gamma_type",
                            "n.B.S", "n.B.D", "n.B.T", 
                            "bias_1", "bias_2", "bias_3", 
                            "bias_PB_S", "bias_PB_D", "bias_PB_T",
                            "bias_NB_S", "bias_NB_D", "bias_NB_T",
                            "MSE_PB_S", "MSE_PB_D", "MSE_PB_T",
                            "MSE_NB_S", "MSE_NB_D", "MSE_NB_T"
)
}

print(Sys.time()-time.start)
write.csv(final.out.mat, paste0(scen.all.ind, "_results_boot_", n.itt, ".csv"))






