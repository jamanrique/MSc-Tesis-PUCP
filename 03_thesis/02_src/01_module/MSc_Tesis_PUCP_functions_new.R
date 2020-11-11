rm(list=ls())

ct <- function(alpha,t){(-log(1-t))^(1/alpha)}
qt <- function(sigma, alpha, t){sigma * ct(alpha,t)}
ptweibull2 <- function(y,qt,alpha,t){1-exp((-ct(alpha,t))*((y/qt)^(alpha)))}
# ptweibull <- function(value, qt, sigma, t) {
  alpha = 1 / (log(sigma + 1))
  ct_2 = ct(alpha,t)
  b = qt / ct_2
  pt = pweibull(value, alpha, b)
  return(pt)
}
qt_2 <- function(df, betas){exp(df %*% betas)}
lancet <- function(){
  temp <- tempfile()
  download.file("http://iinei.inei.gob.pe/iinei/srienaho/descarga/SPSS/447-Modulo552.zip", temp)
  salud <- read_sav(unz(temp, "447-Modulo552/04_C2_CAPITULOS.sav"))
  unlink(temp)
  
  ## pre-procesamiento de datos ##
  salud <- salud[salud$C2P28 != 7,] ## Eliminación de no respuesta
  ## DATA MANAGEMENT
  variables <- c("C2P1","C2P4","C2P7","C2P11","C2P12","C2P13","C2P22A","C2P24","C2P26","REGION","INSTITUCION","C2P28")
  newsalud <- salud[variables]
  newsalud <- data.frame(newsalud[complete.cases(newsalud),])
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
  ## Converting to factor and removing unused tag
  newsalud$C2P7 <- factor(newsalud$C2P7,labels = names(attributes(newsalud$C2P7)$labels)[1:6])
  newsalud$exp <- 2015-newsalud$C2P12
  newsalud <- newsalud %>% dplyr::select(-c(C2P28,C2P12))
  
  ## Converting to factor the remaining categorical variables
  for (j in 1:13) {
    if (j != 6 & j != 13 & j != 11 & j !=12 & j != 3) {
      newsalud[, j] <- factor(newsalud[, j],
                              labels = names(attributes(newsalud[, j])$labels))}}
  
  # levels(newsalud$C2P7)[levels(newsalud$C2P7) == "Locaci?n de servicios (Honorarios profesionales)" | levels(newsalud$C2P7) == "Contrato Administrativo de Servicios (CAS)"] <- "Plazo no fijo"
  # levels(newsalud$C2P7)[levels(newsalud$C2P7) == "Contrato a plazo fijo (sujeto a modalidad)"| levels(newsalud$C2P7) == "Nombrado, permanente" | levels(newsalud$C2P7) == "Plazo indeterminado o indefinido (D.S.728)"] <- "Plazo fijo"
  # 
  # ## C2P23
  # levels(newsalud$C2P23)[levels(newsalud$C2P23) == "Permanente (Tiene trabajo durante todo el a?o de manera continua)?"] <- "Permanente"
  # levels(newsalud$C2P23)[levels(newsalud$C2P23) == "Temporal o estacional (No permanente)?"] <-"Temporal"
  
  ## Re-ordering levels
  newsalud$C2P4 <- relevel(newsalud$C2P4, ref = "Hombre")
  # newsalud$C2P13 <- relevel(newsalud$C2P13, ref = "No")
  # newsalud$C2P23 <- relevel(newsalud$C2P23, ref = "Temporal")
  # newsalud$C2P24 <- relevel(newsalud$C2P24, ref = "No")
  # newsalud$C2P25 <- relevel(newsalud$C2P25, ref = "No")
  # newsalud$C2P26 <- relevel(newsalud$C2P26, ref = "No")
  
  ## Generating separated physicians and nurses data
  
  newsalud <- newsalud %>% dplyr::select(C2P4,li,lf)
  return(newsalud)
}
data_mgmt <- function(data, li, lf){
  BD_l_inf = data[,li]
  BD_l_sup = data[,lf]
  covar = as.data.frame(data[, -c(li, lf)])
  covar = model.matrix(~., covar)
  covar_sin_intercepto = covar[,2:ncol(covar)]
  return(list(lim.inferior = BD_l_inf,
              lim.superior = BD_l_sup,
              matriz.diseno = covar,
              matriz.nointercepto = covar_sin_intercepto))}

ll = function(param,Betas,dbf,t){
  n = length(param) - 1
  for (i in 1:n){Betas[i]=param[i]}
  alpha = param[length(param)]
  qt= qt_2(dbf[['matriz.diseno']], Betas) 
  -sum(log(ptweibull2(dbf[['lim.superior']],qt,alpha,t)-ptweibull2(dbf[['lim.inferior']],qt,alpha,t)))
  # -sum(logDiffExp(ptweibull(dbf[['lim.superior']],qt,alpha,t),ptweibull(dbf[['lim.inferior']],qt,alpha,t)))
}

reg_ces_wei = function(data, li, lf, t, param) {
  dbf <- data_mgmt(data,li,lf)
  Betas <- as.matrix(rep(0,ncol(dbf[['matriz.diseno']])))
  fit_mv <- nlminb(start = param, objective = ll, hessian = T, Betas = Betas, dbf = dbf, t = t)
  # param <- fit_mv$par
  # fit_mv <- optim(par = param, fn = ll, Betas = Betas, dbf = dbf, t = t, hessian = T, method = 'L-BFGS-B',
          # lower = c(rep(-Inf, ncol(dbf[[3]])), 0.10), upper = c(rep(Inf, ncol(dbf[[3]])+1)))
  print(fit_mv$message)
  return(fit_mv)
}
lancet2 <- function(){
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
  
  newsalud_enf <- newsalud_enf %>% dplyr::select(C2P21, C2P27, C2P9, sexo, li,lf)
  return(newsalud_enf)
}

## Tests ##
test <- lancet2()
test$li <- test$li/1000
test$lf <- test$lf/1000
theta_inf <- 5
theta_sup <- 6
t_sim = seq(0.1,0.9,0.1)

gen.cens(family="WEI3",type="interval")
m0 <- gamlss(Surv(li,lf,type="interval2")~.,family = WEI3ic,data = test)
init <- as.vector(c(m0[["mu.coefficients"]],m0[["sigma.coefficients"]]));init

var = list()
for (j in 1:length(t_sim)) {var = append(var,list(reg_ces_wei(data = test,theta_inf,theta_sup,t_sim[j],init)))}

val <- matrix(ncol=length(var[[1]]$par),nrow = 0)
for (l in 1:length(t_sim)) {pba <- var[[l]];val <- rbind(val,pba$par)}
exp(val)

