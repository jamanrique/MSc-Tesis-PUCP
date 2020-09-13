# Justo Manrique - 20091107
# Regresi?n cuant?lica en datos intervalares: un estudio sobre la disparidad de ingresos entre hombres y mujeres en el sistema de salud peruano

#### Carga de librer?as ####
library(tidyverse)
library(haven)
library(dplyr)
library(kyotil)
library(ggpubr)
library(ggthemes)
library(BB)
library(nloptr)
library(optimx)
library(optextras)

#### Segunda parte: creaci?n de la funci?n de regresi?n #####

info <- lancet()

fit_sim_0.5= list()

M = 1000
t_sim = 0.5
b_1 = 0.3
b_2 = 0.6
b_3 = 0.8
sig = 2

X2 = rbeta(10000,2,3)
X3 = rnorm(10000,3,0.5)
X4 = rgamma(10000,2,4)
Qt  = exp(b_1*X2 + b_2*X3 + b_3*X4)

for(j in 1:length(t_sim)){
    Yh = c()
    for(k in 1:M){
        for(i in 1:length(Qt)){Yh[i]= rtweibull(1,Qt[i],sig,t_sim[j])}
        sim = data.frame(
                Li = ifelse(Yh < quantile(Yh,0.25),0,
                        ifelse(Yh <  quantile(Yh,0.5),quantile(Yh,0.25),
                        ifelse(Yh < quantile(Yh,0.75),quantile(Yh,0.50),quantile(Yh,0.70)))),
                Ls = ifelse(Yh < quantile(Yh,0.25),quantile(Yh,0.25),
                        ifelse(Yh <  quantile(Yh,0.5),quantile(Yh,0.5),
                        ifelse(Yh < quantile(Yh,0.75),quantile(Yh,0.75),quantile(Yh,0.99)))),
                X2 =X2,
                X3 = X3,
                X4 = X4)
        fit_sim_0.5  <- append(fit_sim_0.5, list(reg_ces_wei(sim,1,2,t_sim[j])))
        print(paste("Resultado final de la simulación Nª ",k,", para el cuantil ",t_sim[j], ": ",fit_sim_0.5[[k]]$message, sep = ""))
        print(paste("Valores óptimos de la ",k," optimización: ",sep=""))
        print(fit_sim_0.5[[k]]$par)
    }
}

crit=qnorm(0.975)

coverage = list(b_1 = list(),b_2=list(),b_3=list(),sig=list())
sesgo = c()
ecm = c()

for(l in 1:M){
    pba = fit_sim_0.5[[l]]
    
    b_0_test = pba$par[1] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[1])
    b_1_test = pba$par[2] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[2])
    b_2_test = pba$par[3] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[3])
    b_3_test = pba$par[4] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[4])
    sig_test = pba$par[5] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[5])
    
    if (between(b_1,b_1_test[1],b_1_test[2])) {coverage[[1]] = append(coverage[[1]],1)}
    if (between(b_2,b_2_test[1],b_2_test[2])) {coverage[[2]] = append(coverage[[2]],1)}
    if (between(b_3,b_3_test[1],b_3_test[2])) {coverage[[3]] = append(coverage[[3]],1)}
    if (between(sig,sig_test[1],sig_test[2])) {coverage[[4]] = append(coverage[[4]],1)}
    sesgo = append(sesgo,sum(pba$par[2:5] - c(b_1,b_2,b_3,sig)))
    ecm = append(ecm,sum(pba$par[2:5] - c(b_1,b_2,b_3,sig))**2)
}
