## Maestría en Estadística - PUCP ##
## Nombre: Justo A. Manrique Urbina ##
## Modulo MSc_Tesis_PUCP ##

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

rtweibull = function(n, qt, sigma, t) {
  ## Reparametrizaci?n de los datos a la f?rmula b?sica de Weibull en R ##
  alpha = 1 / (log(sigma + 1))
  beta = qt / (-log(1 - t)) ^ (1 / alpha)
  qm = rweibull(n = n, scale = beta, shape = alpha)
  return(qm)
}

dtweibull <- function(value, qt, alpha, t) {
  if (t > 1 | t < 0) {stop("par?metro t solo puede contener valores entre 0 y 1")}
  if (qt < 0) {stop("par?metro qt solo puede ser superior a 0")}
  if (alpha < 0) {stop("par?metro alpha solo puede ser superior a 0")}
  ct = (-log(1 - t)) ^ (1 / alpha)
  b = qt / ct
  a = alpha
  dt = dweibull(x = value, shape = a,scale =  b)
  return(dt)
}

vtweibull <- function(qt, alpha, t) {
  if (t > 1 | t < 0) {stop("par?metro t solo puede contener valores entre 0 y 1")}
  if (qt < 0) {stop("par?metro qt solo puede ser superior a 0")}
  if (alpha < 0) {stop("par?metro alpha solo puede ser superior a 0")}
  ct = (-log(1 - t)) ^ (1 / alpha)
  var = ((qt ^ 2) / (ct) ^ (1 / alpha)) * (gamma(1 + (2 / alpha)) - (gamma(1 + (1 / alpha))) ^ 2)
  return(var)
}

mtweibull <- function(qt, alpha, t) {
  if (t > 1 | t < 0) {stop("par?metro t solo puede contener valores entre 0 y 1")}
  if (qt < 0) {stop("par?metro qt solo puede ser superior a 0")}
  if (alpha < 0) {stop("par?metro alpha solo puede ser superior a 0")}
  ct = (-log(1 - t)) ^ (1 / alpha)
  mt = ((qt) / (ct ^ (1 / alpha))) * gamma(1 + 1 / alpha)
  return(mt)
}

ptweibull <- function(value, qt, sigma, t) {
  # if (t > 1 | t < 0) { stop("par?metro t solo puede contener valores entre 0 y 1") }
  # if (qt <= 0) { stop("par?metro qt solo puede ser superior a 0") }
  # if (sigma <= 0) { stop("par?metro sigma solo puede ser superior a 0") }
  alpha = 1 / (log(sigma + 1))
  ct = (-log(1 - t)) ^ (1 / alpha)
  b = qt / ct
  a = alpha
  pt = pweibull(value, a, b)
  return(pt)
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
  newsalud_enf <- dplyr::select(newsalud_enf, C2P4, li, lf)
  newsalud_enf <- dplyr::mutate(newsalud_enf, pba = if_else(newsalud_enf$C2P4 == "Mujer", 1, 0)) %>% dplyr::select(pba, li, lf)
  return(newsalud_enf)
}

reg_ces_wei = function(data, li, lf, t) {
  data_l_inf = data[,li]
  data_l_sup = data[, lf]
  covar = as.matrix(data[, - c(li, lf)])
  n = ncol(covar)
  d = nrow(covar)
  
  ll = function(param, t) {
    b_0 = 0
    b = as.matrix(rep(0, ncol(covar)))
    b_0 = param[1]        
    for (i in 1:n) {b[i] = param[i + 1]}
    sigma = param[length(param)]
    ct = -(log(1-t))
    qt = exp(b_0 + covar %*% b)
    s_ll = 0
    for (j in 1:d) {
      qt_j = qt[j,]
      s_ll = s_ll - logDiffExp(ptweibull(data_l_sup[j], qt_j, sigma, t), ptweibull(data_l_inf[j], qt_j, sigma, t))
    }
    return(s_ll)
  }
  inicial = c(as.vector(coef(lm(data_l_inf~covar,data=data.frame(cbind(data_l_inf,covar))))),1)
  fit_mv = nlminb(inicial,ll,t = t,lower = c(rep(0,ncol(covar)+1,0.1)),upper =c(rep(Inf,ncol(covar)+2)))
  inicial = fit_mv$par
  print(paste("Resultado de la 1ª Optimización: ",fit_mv$message))
  fit_mv =   optim(par = inicial, fn = ll,t = t, upper = c(rep(Inf, ncol(covar) + 2)),lower = c(rep(0, ncol(covar) + 1), 0.01),hessian = T)
  
  return(fit_mv)
}

