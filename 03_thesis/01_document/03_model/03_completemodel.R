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

fit_sim= list(
    t_0.25 = list(),
    t_0.50 = list(),
    t_0.75 = list()
)


t_sim = seq(0.25,0.75,0.25)
b_1 = 0.3
b_2 = 0.6
b_3 = 0.8
sig = 2

X2 = rbeta(10000,2,3)
X3 = rnorm(10000,2,0.25)
X4 = rgamma(10000,2,25)
Qt  = exp(b_1*X2 + b_2*X3 + b_3*X4)

for(j in 1:length(t_sim)){
    Yh = c()
    M = 2
    while (M != 0) {
        try({
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
            fit_sim[[j]]  <- append(fit_sim[[j]], list(reg_ces_wei(sim,1,2,t_sim[j])))
            print(paste("Resultado final de la simulación Nª ",abs(2-M+1),", para el cuantil ",t_sim[j], ": CONVIRGIÓ", sep = ""))
            M <- M - 1
            }
        )
    M=2
    }
}

crit=qnorm(0.975)

coverage = list(
    t_0.25 = list( b_1 = list(),b_2=list(),b_3=list(),sig=list()),
    t_0.50 = list( b_1 = list(),b_2=list(),b_3=list(),sig=list()),
    t_0.75 = list( b_1 = list(),b_2=list(),b_3=list(),sig=list())
    )

sesgo = list(
    t_0.25 = list(),
    t_0.50 = list(),
    t_0.75 = list()
)
ecm = list(
    t_0.25 = list(),
    t_0.50 = list(),
    t_0.75 = list()
)

for(j in 1:length(t_sim)){
    for(l in 1:M){
        pba = fit_sim[[j]][[l]]
        b_0_test = pba$par[1] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[1])
        b_1_test = pba$par[2] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[2])
        b_2_test = pba$par[3] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[3])
        b_3_test = pba$par[4] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[4])
        sig_test = pba$par[5] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[5])
    
        if (between(b_1,b_1_test[1],b_1_test[2])) {coverage[[j]][[1]] = append(coverage[[j]][[1]],1)}
        if (between(b_2,b_2_test[1],b_2_test[2])) {coverage[[j]][[2]] = append(coverage[[j]][[2]],1)}
        if (between(b_3,b_3_test[1],b_3_test[2])) {coverage[[j]][[3]] = append(coverage[[j]][[3]],1)}
        if (between(sig,sig_test[1],sig_test[2])) {coverage[[j]][[4]] = append(coverage[[j]][[4]],1)}
        sesgo[[j]] = append(sesgo[[j]],sum(pba$par[2:5] - c(b_1,b_2,b_3,sig)))
        ecm[[j]] = append(ecm[[j]],sum(pba$par[2:5] - c(b_1,b_2,b_3,sig))**2)
}
}
