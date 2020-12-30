rm(list=ls())
library(gamlss)
library(gamlss.cens)
library(haven)
library(tidyverse)
library(BB)
library(matrixStats)
library(gridExtra)
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
  fit_mv = optim(par = param,fn = likelihood,hessian = T,betas=Bs,df=df,tau=tau)
  print(fit_mv$convergence)
  return(fit_mv)
}
lancet <- function(){
  temp <- tempfile()
  download.file("http://iinei.inei.gob.pe/iinei/srienaho/descarga/SPSS/447-Modulo552.zip", temp)
  salud <- read_sav(unz(temp, "447-Modulo552/04_C2_CAPITULOS.sav"))
  unlink(temp)
  
  ## Eliminating missing values
  salud <- salud[salud$C2P28 != 7,]
  ## DATA MANAGEMENT
  myvarsalud <- c("C2P21", "C2P27", "C2P9", "REGION", "INSTITUCION", "C2P1", "C2P4","C2P7", "C2P11", "C2P13", "C2P23", "C2P24", "C2P25", "C2P26", "C2P28")
  
  newsalud <- salud[myvarsalud]
  newsalud <- data.frame(newsalud)
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
  
  ## Deleting "Other" in type of contract
  newsalud <- newsalud[newsalud$C2P7 != 6,]
  ## Deleting 8 NA's in dependents
  newsalud <- newsalud[!is.na(newsalud$C2P9),]
  ## Converting to factor and removing unused tag
  newsalud$C2P28 <- factor(newsalud$C2P28,labels = names(attributes(newsalud$C2P28)$labels)[1:6])
  newsalud$C2P7 <- factor(newsalud$C2P7,labels = names(attributes(newsalud$C2P7)$labels)[1:5])
  
  ## Converting to factor the remaining categorical variables
  for (j in 4:15) {
    if (j != 8 & j != 15) {
      newsalud[, j] <- factor(newsalud[, j],
                              labels = names(attributes(newsalud[, j])$labels))}}
  
  ## REGION
  levels(newsalud$REGION)[levels(newsalud$REGION) == "Sierra"] <- "Nuevo"
  levels(newsalud$REGION)[levels(newsalud$REGION) == "Selva"] <- "Sierra"
  levels(newsalud$REGION)[levels(newsalud$REGION) == "Nuevo"] <- "Selva"
  
  ## C2P7
  levels(newsalud$C2P7)[levels(newsalud$C2P7) == "Locaci?n de servicios (Honorarios profesionales)" | levels(newsalud$C2P7) == "Contrato Administrativo de Servicios (CAS)"] <- "Plazo no fijo"
  levels(newsalud$C2P7)[levels(newsalud$C2P7) == "Contrato a plazo fijo (sujeto a modalidad)"| levels(newsalud$C2P7) == "Nombrado, permanente" | levels(newsalud$C2P7) == "Plazo indeterminado o indefinido (D.S.728)"] <- "Plazo fijo"
  
  ## C2P23
  levels(newsalud$C2P23)[levels(newsalud$C2P23) == "Permanente (Tiene trabajo durante todo el a?o de manera continua)?"] <- "Permanente"
  levels(newsalud$C2P23)[levels(newsalud$C2P23) == "Temporal o estacional (No permanente)?"] <-"Temporal"
  
  ## Re-ordering levels
  newsalud$C2P4 <- relevel(newsalud$C2P4, ref = "Mujer")
  newsalud$C2P13 <- relevel(newsalud$C2P13, ref = "No")
  newsalud$C2P23 <- relevel(newsalud$C2P23, ref = "Temporal")
  newsalud$C2P24 <- relevel(newsalud$C2P24, ref = "No")
  newsalud$C2P25 <- relevel(newsalud$C2P25, ref = "No")
  newsalud$C2P26 <- relevel(newsalud$C2P26, ref = "No")
  
  ## Generating separated physicians and nurses data
  
  newsalud <- subset(newsalud, select = -c(C2P28))
  newsalud_med <- newsalud[newsalud$C2P1 == "M?dico",]
  newsalud_med <- subset(newsalud_med, select = -C2P1)
  newsalud_enf <- newsalud[newsalud$C2P1 != "M?dico",]
  newsalud_enf <- subset(newsalud_enf, select = -C2P1)
  
  ## Scaling
  newsalud_med$li <- newsalud_med$li;
  newsalud_med$lf <- newsalud_med$lf
  
  newsalud_enf$li <- newsalud_enf$li
  newsalud_enf$lf <- newsalud_enf$lf
  
  newsalud_enf <- subset(newsalud_enf, (newsalud_enf$lf == Inf) == FALSE)
  newsalud_enf <- dplyr::mutate(newsalud_enf, 
                                sexo = if_else(newsalud_enf$C2P4 == "Mujer", 1, 0)
  )
  
  newsalud_enf <- newsalud_enf %>% dplyr::select(C2P21,C2P27,C2P9,sexo, li,lf)
  newsalud_enf$sexo <- as.factor(newsalud_enf$sexo)
  return(newsalud_enf)
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

final_database

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
round(final_database,4)

#### Datos reales ####

real_data = lancet()

tau_seq_sim = seq(0.10,0.9,0.1)
m1 = gamlss(Surv(li,lf,type="interval2")~.,family = WEI3ic,data = real_data)
init_real = as.vector(c(coef(m1),m1$sigma.coefficients));init_real

var = list()
for (j in 1:length(tau_seq_sim)) {
  var = append(var,list(reg_Wr(data = real_data,li = 5, lf= 6,tau_seq_sim[j],param = init_real)))
}  

val <- matrix(ncol=length(var[[1]]$par),nrow = 0)
for (l in 1:length(tau_seq_sim)) {
  pba <- var[[l]]
  val <- rbind(val,pba$par)
}
val
pred <- as.data.frame(rbind(var[[1]]$par - qnorm(0.975) * sqrt(diag(solve(var[[1]]$hessian))),
                            var[[1]]$par + qnorm(0.975) * sqrt(diag(solve(var[[1]]$hessian)))))

val
ggplot_aa <- list()
for (k in 1:6) {
  ar <- matrix(nrow=0,ncol=2)
  for (t in 1:length(tau_seq_sim)) {
   ar <- rbind(ar,cbind(var[[t]]$par[k]  + qnorm(0.975) * sqrt(diag(solve(var[[t]]$hessian))[k]),
         var[[t]]$par[k]  - qnorm(0.975) * sqrt(diag(solve(var[[t]]$hessian))[k])
   ))}
  ar <- as.data.frame(ar)

  p = ggplot(data=as.data.frame(val),
              aes_string(x=tau_seq_sim,y=val[,k])) +
              geom_point() +
              geom_line() +
              geom_ribbon(data = ar,
                          aes_string(ymin=ar[,1],ymax=ar[,2]),
                          alpha=0.4)
  ggplot_aa[[k]] = p
  }

plot_grid_alpha <- do.call("grid.arrange",ggplot_aa)
