rm(list=ls(all=TRUE))
library(foreign)
library(nnet)
library(gmodels)
library(VGAM)
library(rstan)
library(haven)    
library(lme4)
library(tidyverse)
salud.personal = read_sav("http://portal.susalud.gob.pe/wp-content/uploads/archivo/base-de-datos/2015/CUESTIONARIO%2002%20-%20CAPITULOS.sav")
salud.medicos_2  =  salud.personal[salud.personal$C2P1==1,]
salud.medicos$C2P4 = factor(salud.medicos$C2P4,levels=attr(salud.medicos_2$C2P4,"labels"))
salud.medicos$C2P13 = factor(salud.medicos$C2P13,levels=attr(salud.medicos_2$C2P13,"labels"))
salud.medicos$C2P25 = factor(salud.medicos$C2P25,levels=attr(salud.medicos_2$C2P25,"labels"))
salud.medicos$C2P26 = factor(salud.medicos$C2P26,levels=attr(salud.medicos_2$C2P26,"labels"))

model.nom = multinom(INSTITUCION ~ C2P4+C2P9+C2P2EDAD+C2P21+C2P27+C2P13+C2P24+C2P25+C2P26,salud.medicos)
summary(model.nom)

exp(coef(model.nom))

salud.medicos_2$INSTITUCION
salud.medicos_2$C2P13
salud.medicos_2$C2P25

salud.medicos     = salud.medicos[salud.medicos$C2P28 < 7,]
salud.medicos$I   = ifelse(salud.medicos$C2P28==1,1,0)
salud.medicos$II  = ifelse(salud.medicos$C2P28==2,1,0)
salud.medicos$III = ifelse(salud.medicos$C2P28==3,1,0)
salud.medicos$IV  = ifelse(salud.medicos$C2P28==4,1,0)
salud.medicos$V   = ifelse(salud.medicos$C2P28==5,1,0)
salud.medicos$VI  = ifelse(salud.medicos$C2P28==6,1,0)

model.1  = vglm(cbind(I,II,III,IV,V,VI) ~ factor(INSTITUCION) + C2P4 + C2P9 + C2P2EDAD + C2P21+ C2P27 + C2P13 + C2P24+ C2P25 + C2P26,
                family = acat(parallel = TRUE),
                data=salud.medicos)
summary(model.1)
exp(coef(model.1))

model.2  = vglm(cbind(I,II,III,IV,V,VI) ~ factor(INSTITUCION) + C2P4 + C2P9 + C2P2EDAD + C2P21+ C2P27 + C2P13 + C2P24+ C2P25 + C2P26,
                family = cumulative(parallel = TRUE),
                data=salud.medicos)
summary(model.2)
round(exp(coef(model.2)),4)

model.3  = vglm(cbind(I,II,III,IV,V,VI) ~ factor(INSTITUCION) + C2P4 + C2P9 + C2P2EDAD + C2P21+ C2P27 + C2P13 + C2P24+ C2P25 + C2P26,
                family = sratio(reverse=TRUE, parallel = TRUE),
                data=salud.medicos)
summary(model.3)
exp(coef(model.3))

lrtest_vglm(model.1,model.2)
lrtest_vglm(model.2,model.3)


true_values = c()
for (m in 1:1000) {
  N     = 20
  q     = 100
  x1     = rnorm(N*q,0,1)
  beta0 = -6
  beta1 = 0.5
  b0    = rnorm(N,0,4)    # Efecto aleatorio por familia
  b0    = rep(b0,times=q)
  p     = exp(beta0 + b0 +  beta1*x1)*(exp(beta0+b0 + beta1*x1) + 1)^(-1)
  y     = rbinom(N*q,1,prob = p)
  datos = data.frame(id=1:N*q,y=y,x1=x1,
                     familia=rep(1:N,5))
  
  model.glm = glmer(y ~ x1 + (1| id ),family = binomial(),data = datos)
  conf = confint(model.glm, method = "boot",nsim = 10)
  true_values=append(true_values,between(beta1,conf[3,1],conf[3,2]))
}

sum(true_values)/length(true_values)


