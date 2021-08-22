library(gamlss)
library(gamlss.cens)
library(haven)
library(tidyverse)
library(BB)
library(matrixStats)
library(gridExtra)
library(numDeriv)
library(pracma)
library(ucminf)
library(nloptr)
library(ggthemes)

gen.cens(family="WEI3",type="interval")

Ct = function(alpha, tau){
  (-log(1-tau))^(1/alpha)
}
F_Wr = function(Y,Qt,alpha,tau){
  B = Qt/Ct(alpha,tau)
  pweibull(Y,shape=alpha,scale=B)
}
Qt_b = function(betas, df){
  Qt_c = exp(as.matrix(df[['matriz.diseno']]) %*% betas)
  return(Qt_c)
}
Qt_a = function(betas, df){
  Qt_c = exp(as.matrix(df) %*% betas)
  return(Qt_c)
}
Rand_Wr = function(n,Qt,alpha,tau){
  B = Qt/Ct(alpha,tau)
  rweibull(n,shape = alpha,scale = B)
}
Simulation = function(n){
  X1 = rnorm(n,2,0.25)
  X2 = rbeta(n,2,3)
  X3 = rgamma(n,2,20)
  df = data.frame(X1 = X1, X2 = X2, X3 = X3)
  return(df)
}
DF_Simulation = function(df,betas,alpha,tau){
  M = dim(df)[1]
  design_matrix = model.matrix(~ . ,df)
  Qt_i = Qt_a(betas,design_matrix)
  Y = c()
  for (j in 1:M) {
    Y=rbind(Y,Rand_Wr(1,Qt_i[j],alpha,tau))
  }
  min_Y = floor(min(Y))
  if (min_Y==0) {
    min_Y <- 0.01
  }
  Q8_Y = round(quantile(Y,0.8),2)
  interval = round((Q8_Y-min_Y) / 6,2)
  seq_interv = c(seq(min_Y,Q8_Y,interval),Inf)
  Ls = c()
  for (u in 1:length(Y)) {
    for (n in 1:length(seq_interv)) {
      if (Y[u] <seq_interv[n]) {
        Ls[u] = seq_interv[n]
        break
      }}}
  Li = c()
  for (p in 1:length(Y)) {
    for (w in 1:length(seq_interv)) {
      if (Y[p] > rev(seq_interv)[w]) {
        Li[p] = rev(seq_interv)[w]
        break
      }}}
  F_df = cbind(data.frame(Li = round(Li,2), Ls = round(Ls,2)),df)
  return(F_df)
}
data_mgmt = function(data,li,lf){
  df_inf = data[,li]
  df_sup = data[,lf]
  covar = as.data.frame(data[,-c(li,lf)])
  covar = model.matrix(~. , covar)
  lista = list(lim.inferior = df_inf, lim.superior = df_sup, matriz.diseno = covar)
  return(lista)
}
likelihood = function(param,betas,df,tau){
  n = length(param)-1
  for(i in 1:n){betas[i]=param[i]}
  alpha = param[length(param)]
  Qt_i = Qt_b(betas,df)
  -sum(log(F_Wr(df[['lim.superior']],Qt_i,alpha,tau)-F_Wr(df[['lim.inferior']],Qt_i,alpha,tau)))
}
reg_Wr = function(data,li,lf,tau,param){
  df = data_mgmt(data,li,lf)
  Bs = as.matrix(rep(0,ncol(df[['matriz.diseno']])))
  fit_mv = nloptr(x0 = param,eval_f = likelihood, betas=Bs, df=df, tau=tau,opts = list("algorithm"=c("NLOPT_LN_SBPLX"),maxeval = 400),lb = c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf),ub=c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf))
  print(fit_mv$message)
  fit_mv$hessian = pracma::hessian(x0 = fit_mv$solution,f = likelihood,betas=Bs,df=df,tau=tau)
  return(fit_mv)
}

lancet <- function(int){
  temp <- tempfile()
  download.file("http://iinei.inei.gob.pe/iinei/srienaho/descarga/SPSS/447-Modulo552.zip", temp)
  salud <- read_sav(unz(temp, "447-Modulo552/04_C2_CAPITULOS.sav"))
  unlink(temp)
  
  ## Manteniendo respuestas válidas
  salud <- salud[salud$C2P28 != 7,]
  
  salud$C2P26
  
  ### Pre-procesamiento ###
  variables <- c("INSTITUCION","C2P4","C2P13","C2P21","C2P24","C2P26","C2P28","C2P1","C2P27")
  
  newsalud <- salud[variables]
  
  ## Disaggregating salary ordinal variable into its limits
  newsalud$li <- NULL
  newsalud$ls <- NULL
  newsalud$li <- ifelse(newsalud$C2P28 == 1, 850,
                        ifelse(newsalud$C2P28 == 2, 1000,
                               ifelse(newsalud$C2P28 == 3, 2001,
                                      ifelse(newsalud$C2P28 == 4, 3001,
                                             ifelse(newsalud$C2P28 == 5, 4001, 5001)))))
  newsalud$lf <- ifelse(newsalud$C2P28 == 1, 999,
                        ifelse(newsalud$C2P28 == 2, 2000,
                               ifelse(newsalud$C2P28 == 3, 3000,
                                      ifelse(newsalud$C2P28 == 4, 4000,
                                             ifelse(newsalud$C2P28 == 5, 5000, Inf)))))
  
  newsalud <- newsalud %>% select(-C2P28)
  
  newsalud$INSTITUCION <- factor(newsalud$INSTITUCION,labels=names(attr(newsalud$INSTITUCION,"labels")))
  newsalud$C2P4 <- factor(newsalud$C2P4,labels = names(attr(newsalud$C2P4,"labels")))
  newsalud$C2P13 <- factor(newsalud$C2P13, labels=names(attr(newsalud$C2P13,"labels")))
  newsalud$C2P21 <- as.integer(newsalud$C2P21)
  newsalud$C2P24 <- factor(newsalud$C2P24, labels = names(attr(newsalud$C2P24,"labels")))
  newsalud$C2P26 <- factor(newsalud$C2P26, labels = names(attr(newsalud$C2P26,"labels")))
  newsalud$C2P1 <- factor(newsalud$C2P1, labels = names(attr(newsalud$C2P1,"labels")))
  newsalud$C2P27 <- as.integer(newsalud$C2P27)
  
  ## Generating separated physicians and nurses data
  
  newsalud_enf <- newsalud[newsalud$C2P1 =="Enfermero/a" ,]
  newsalud_enf <- subset(newsalud_enf, select = -C2P1)
  newsalud_enf <- as.data.frame(newsalud_enf)
  
  newsalud_med <- newsalud[newsalud$C2P1 == "Médico",]
  newsalud_med <- subset(newsalud_med, select = -C2P1)
  newsalud_med <- as.data.frame(newsalud_med)
  if (int == 1) {
    return(newsalud_enf)  
  }
  if (int ==2) {
    return(newsalud_med)
  }
  if (int ==3) {
    newsalud_total <- as.data.frame(subset(newsalud,select = -C2P1))
    return(newsalud_total)
  }
}

#### Simulación ####

L = 5000
n = c(100,500,1000)
betas_sim = c(7,0.3,0.84,2.5)
alpha_sim = 2
tau_sim = seq(0.1,0.9,0.1)

sim_list= list(n100 = list(), n500 = list(), n1000 = list())

for (j in 1:length(n)) {
  for (k in 1:length(tau_sim)) {
    for (p in 1:L) {
      sim = Simulation(n[j])
      df_sim = DF_Simulation(sim,betas_sim,alpha_sim,tau_sim[k])
      m0 = gamlss(Surv(Li,Ls,type="interval2")~.,family = WEI3ic,data = na.omit(df_sim))
      init = as.vector(c(coef(m0),m0$sigma.coefficients))
      sim_list[[j]] = append(sim_list[[j]],list(reg_Wr(data = df_sim,li = 1,lf = 2,tau = tau_sim[k],param = init)))
    }
  }
}

se_calc <- function(lista) {
  ic_contain <- matrix(nrow=0,ncol=5)
  for (k in 1:length(lista)) {
    little_t <- lista[[k]]
    ic <- rbind( little_t$par - qnorm(0.975) * sqrt(diag(solve(little_t$hessian))),
                 little_t$par + qnorm(0.975) * sqrt(diag(solve(little_t$hessian))))
    ic_contain <- rbind(ic_contain,cbind(between(betas_sim[1],left = ic[1,1],right = ic[2,1]),
                        between(betas_sim[2],left = ic[1,2],right = ic[2,2]),
                        between(betas_sim[3],left = ic[1,3],right = ic[2,3]),
                        between(betas_sim[4],left = ic[1,4],right = ic[2,4]),
                        between(alpha_sim,left = ic[1,5],right = ic[2,5])))
  }
  final_df <- as.matrix(t(colSums(ic_contain)/length(lista)))
  return(final_df)
}

final_database <- matrix(nrow=0,ncol=5)
seq_t <- seq(1,8,1);seq_t

for (j in 1:length(n)) {
  df <- sim_list[[j]]
  actual_t <- df[1:5000*seq_t[1]]
  vf_df <- se_calc(actual_t)
  rownames(vf_df) <- paste("n_",n[j],sep="")
  final_database <- rbind(final_database,vf_df)
  
  for (l in 1:length(seq_t)) {
   actual_t <-  df[(5000*seq_t[l]+1):(5000*(seq_t[l]+1))]
   vf_df <- se_calc(actual_t)
   rownames(vf_df) <- paste("n_",n[j],sep="")
   final_database <- rbind(final_database,vf_df)
  }
}

final_database <- as.data.frame(final_database)

final_database_vf <- cbind(
  Cuantil = rep(c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'),3),
  B1 = c(rep(100,9),rep(500,9),rep(1000,9)), 
  final_database)

v1 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V1) +
  geom_boxplot(shape = "circle", fill = "#B22222") +
  facet_wrap(vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Cobertura",title = "Cobertura para \U003B2_0") +
  scale_y_continuous(limits=c(0.93,0.96)) + 
  geom_hline(yintercept=0.95,linetype="dashed",color="red")+
  theme_minimal() +
  theme(legend.position="bottom")

v2 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V2) +
  geom_boxplot(shape = "circle", fill = "#B22222") +
  facet_wrap(vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Cobertura",title = "Cobertura para \U003B2_1") +
  scale_y_continuous(limits=c(0.93,0.96)) + 
  geom_hline(yintercept=0.95,linetype="dashed",color="red")+
  theme_minimal() +
  theme(legend.position="bottom")

v3 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V3) +
  geom_boxplot(shape = "circle", fill = "#B22222") +
  facet_wrap(vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Cobertura",title = "Cobertura para \U003B2_2") +
  scale_y_continuous(limits=c(0.93,0.96)) + 
  geom_hline(yintercept=0.95,linetype="dashed",color="red")+
  theme_minimal() +
  theme(legend.position="bottom")

v4 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V4) +
  geom_boxplot(shape = "circle", fill = "#B22222") +
  facet_wrap(vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Cobertura",title = "Cobertura para \U003B2_3") +
  scale_y_continuous(limits=c(0.93,0.96)) + 
  geom_hline(yintercept=0.95,linetype="dashed",color="red")+
  theme_minimal() +
  theme(legend.position="bottom")

v5 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V5) +
  geom_boxplot(shape = "circle", fill = "#B22222") +
  facet_wrap(vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Cobertura",title = "Cobertura para \U003B1") +
  scale_y_continuous(limits=c(0.93,0.96)) + 
  geom_hline(yintercept=0.95,linetype="dashed",color="red")+
  theme_minimal() +
  theme(legend.position="bottom")

ggpubr::ggarrange(v1,v2,v3,v4,v5)

ecm_calc <- function(lista) {
  ic_contain <- matrix(nrow=0,ncol=5)
  for (k in 1:length(lista)) {
    little_t <- lista[[k]]
    ecm <- (little_t$par - c(betas_sim,alpha_sim))**2
    ic_contain <- rbind(ic_contain,ecm)
  }
  final_df <- as.matrix(t(colSums(ic_contain)/length(lista)))
  return(final_df)
}

  final_database <- matrix(nrow=0,ncol=5)

seq_t <- seq(1,8,1);seq_t

for (j in 1:length(n)) {
  df <- sim_list[[j]]
  actual_t <- df[1:5000*seq_t[1]]
  vf_df <- ecm_calc(actual_t)
  rownames(vf_df) <- paste("n_",n[j],sep="")
  final_database <- rbind(final_database,vf_df)
  
  for (l in 1:length(seq_t)) {
    actual_t <-  df[(5000*seq_t[l]+1):(5000*(seq_t[l]+1))]
    vf_df <- ecm_calc(actual_t)
    rownames(vf_df) <- paste("n_",n[j],sep="")
    final_database <- rbind(final_database,vf_df)
  }
}
round(final_database,4)

final_database <- as.data.frame(final_database)

final_database_vf <- cbind(
  Cuantil = rep(c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'),3),
  B1 = c(rep(100,9),rep(500,9),rep(1000,9)), 
  final_database)

v1 <- ggplot(data = final_database_vf) +
  geom_line(aes(x=B1,y=V1,color=Cuantil,group=Cuantil),size=1.2) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="ECM",title = "ECM para \U003B2_0") +
  scale_x_continuous(breaks=c(100,500,1000)) +
  theme_minimal() +
  theme(legend.position="bottom")


v2 <- ggplot(data = final_database_vf) +
  geom_line(aes(x=B1,y=V2,color=Cuantil,group=Cuantil),size=1.2) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="ECM",title = "ECM para \U003B2_1") +
  scale_x_continuous(breaks=c(100,500,1000)) + 
  theme_minimal() +
  theme(legend.position="bottom")

v3 <- ggplot(data = final_database_vf) +
  geom_line(aes(x=B1,y=V3,color=Cuantil,group=Cuantil),size=1.2) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="ECM",title = "ECM para \U003B2_2") +
  scale_x_continuous(breaks=c(100,500,1000)) + 
  theme_minimal() +
  theme(legend.position="bottom")

v4 <- ggplot(data = final_database_vf) +
  geom_line(aes(x=B1,y=V4,color=Cuantil,group=Cuantil),size=1.2) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="ECM",title = "ECM para \U003B2_3") +
  scale_x_continuous(breaks=c(100,500,1000)) + 
  theme_minimal() +
  theme(legend.position="bottom")

v5 <- ggplot(data = final_database_vf) +
  geom_line(aes(x=B1,y=V5,color=Cuantil,group=Cuantil),size=1.2) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="ECM",title = "ECM para \U003B1") +
  scale_x_continuous(breaks=c(100,500,1000)) + 
  theme_minimal() +
  theme(legend.position="bottom")

ggpubr::ggarrange(v1,v2,v3,v4,v5)

ses_calc <- function(lista) {
  ic_contain <- matrix(nrow=0,ncol=5)
  for (k in 1:length(lista)) {
    little_t <- lista[[k]]
    ecm <- (little_t$par - c(betas_sim,alpha_sim))
    ic_contain <- rbind(ic_contain,ecm)
  }
  final_df <- as.matrix(t(colSums(ic_contain)/length(lista)))
  return(final_df)
}

final_database <- matrix(nrow=0,ncol=5)

seq_t <- seq(1,8,1);seq_t

for (j in 1:length(n)) {
  df <- sim_list[[j]]
  actual_t <- df[1:5000*seq_t[1]]
  vf_df <- ses_calc(actual_t)
  rownames(vf_df) <- paste("n_",n[j],sep="")
  final_database <- rbind(final_database,vf_df)
  
  for (l in 1:length(seq_t)) {
    actual_t <-  df[(5000*seq_t[l]+1):(5000*(seq_t[l]+1))]
    vf_df <- ses_calc(actual_t)
    rownames(vf_df) <- paste("n_",n[j],sep="")
    final_database <- rbind(final_database,vf_df)
  }
}

final_database <- as.data.frame(final_database)

final_database_vf <- cbind(
  Cuantil = rep(c('0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'),3),
  B1 = c(rep(100,9),rep(500,9),rep(1000,9)), 
  final_database)

v1 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V1) +
  geom_boxplot(shape = "circle", fill = "#EF562D") +
  facet_grid(vars(), vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Sesgo Relativo",title = "Sesgo Relativo para \U003B2_0") +
  geom_hline(yintercept=0.0,linetype="dashed",color="black") +
  scale_y_continuous(limits=c(-0.025,0.025)) +
  theme_minimal() +
  theme(legend.position="bottom")

v2 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V1) +
  geom_boxplot(shape = "circle", fill = "#EF562D") +
  facet_grid(vars(), vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Sesgo Relativo",title = "Sesgo Relativo para \U003B2_1") +
  geom_hline(yintercept=0.0,linetype="dashed",color="black") +
  scale_y_continuous(limits=c(-0.025,0.025)) +
  theme_minimal() +
  theme(legend.position="bottom")

v3 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V1) +
  geom_boxplot(shape = "circle", fill = "#EF562D") +
  facet_grid(vars(), vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Sesgo Relativo",title = "Sesgo Relativo para \U003B2_2") +
  geom_hline(yintercept=0.0,linetype="dashed",color="black") +
  scale_y_continuous(limits=c(-0.025,0.025)) +
  theme_minimal() +
  theme(legend.position="bottom")

v4 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V1) +
  geom_boxplot(shape = "circle", fill = "#EF562D") +
  facet_grid(vars(), vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Sesgo Relativo",title = "Sesgo Relativo para \U003B2_3") +
  geom_hline(yintercept=0.0,linetype="dashed",color="black") +
  scale_y_continuous(limits=c(-0.025,0.025)) +
  theme_minimal() +
  theme(legend.position="bottom")

v5 <- ggplot(data = final_database_vf) +
  aes(x = "", y = V1) +
  geom_boxplot(shape = "circle", fill = "#EF562D") +
  facet_grid(vars(), vars(B1)) +
  scale_color_brewer(palette = "YlOrRd") + 
  labs(x="Tamaño de Muestra", y="Sesgo Relativo",title = "Sesgo Relativo para \U003B1") +
  geom_hline(yintercept=0.0,linetype="dashed",color="black") +
  theme_minimal() +
  theme(legend.position="bottom")

ggpubr::ggarrange(v1,v2,v3,v4,v5)

#### Datos reales ####

real_data_enf = lancet(3)
real_data_med = real_data_enf

library(skimr)

enf <- real_data_med %>% mutate(factor = case_when(lf == 999 ~ 1, lf == 2000 ~ 2, lf == 3000 ~ 3, lf == 4000 ~ 4,lf == 5000 ~ 5, lf == Inf ~ 6))
enf$factor <- factor(x = enf$factor,labels = c("[850-999]","[1000-2000]","[2001-3000]","[3001-4000]","[4001-5000]","[5001-Inf.]"))

ib <- enf %>% group_by(factor) %>% skim()

enf %>% group_by(factor,INSTITUCION) %>% summarise( n())



tau_seq_sim = seq(0.1,0.90,0.1)
#m1 = gamlss(Surv(li,lf,type="interval2")~.,family = WEI3ic,data = real_data_enf)
m2 = gamlss(Surv(li,lf,type="interval2")~.,family = WEI3ic,data = real_data_med)
# init_real_enf = as.vector(c(coef(m1),m1$sigma.coefficients))
init_real_med = as.vector(c(coef(m2),exp(m2$sigma.coefficients)))
abc <- confint(m2)
#### Regresión cuantílica para médicxs ####

var_med = list()
for (j in 1:length(tau_seq_sim)) {
  var_med = append(var_med,list(reg_Wr(data = real_data_med,li = 8, lf = 9,tau_seq_sim[j],param = init_real_med)))
}

col_df <- c(colnames(model.matrix( ~ . ,subset.data.frame(real_data_med,select = -c(8,9)))),"\U03B1")
val_med <- matrix(ncol=8,nrow = 0)

for (l in 1:length(tau_seq_sim)) {
  pba <- var_med[[l]]
  val_med <- rbind(val_med,
                   cbind(pba$solution,
                         col_df,
                         tau_seq_sim[l],
                         init_real_med,
                         pba$solution + qnorm(0.975) * sqrt(diag(solve(pba$hessian))),
                         pba$solution - qnorm(0.975) * sqrt(diag(solve(pba$hessian))),
                         c(abc[,1],exp(m2$sigma.coefficients)),
                         c(abc[,2],exp(m2$sigma.coefficients))
                         ))
}
val_med <- as.data.frame(val_med)
val_med$V1 <- exp(as.numeric(val_med$V1))
val_med$V3 <- as.numeric(val_med$V3)
val_med$V5 <- exp(as.numeric(val_med$V5))
val_med$V6 <- exp(as.numeric(val_med$V6))
val_med$V7 <- exp(as.numeric(val_med$V7))
val_med$V8 <- exp(as.numeric(val_med$V8))
val_med$init_real_med <- as.numeric(val_med$init_real_med)


ggplot(val_med) +
  aes(x = V3, y = V1) +
  geom_ribbon(mapping = aes(ymin = V6,ymax = V5),fill = "#B2B9FF") +
  geom_line(size = 0.5, colour = "#1B00FF") +
  geom_ribbon(mapping = aes(ymin = V7,ymax = V8),fill = "#cccccc",alpha=0.5) +
  geom_line(mapping = aes(y=exp(init_real_med)),size=0.5) + 
  ggthemes::theme_pander() +
  facet_wrap(vars(col_df), scales = "free")
