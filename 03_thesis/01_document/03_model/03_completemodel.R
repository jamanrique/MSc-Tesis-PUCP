# Justo Manrique - 20091107
# Regresi?n cuant?lica en datos intervalares: un estudio sobre la disparidad de ingresos entre hombres y mujeres en el sistema de salud peruano

rm(list = ls())
#### Carga de librer?as ####
library(tidyverse)
library(haven)
library(dplyr)
library(kyotil)
library(ggpubr)
library(ggthemes)

#### Creaci?n de la funci?n de distribuci?n para n?meros aleatorios ####
rtweibull = function(n, qt, sigma, t) {
    ## Reparametrizaci?n de los datos a la f?rmula b?sica de Weibull en R ##
    alpha = 1 / (log(sigma + 1))
    beta = qt / (-log(1 - t)) ^ (1 / alpha)
    qm = rweibull(n = n, scale = beta, shape = alpha)
    return(qm)
}

### Funci?n de log-verosimilitud

ll = function(param, x, t) {
    qt <- param[1]
    sigma <- param[2]
    ct <- -(log(1 - t))
    a <- 1 / (log(sigma + 1))

    - sum(log(a) + log(ct) - a * log(qt) + (a - 1) * log(x) - ct * (x / qt) ** a)
}

### Optimizaci?n de la funci?n ###

muestra <- rtweibull(10000, 4, 0.5, 0.5)
hist(muestra)

inicial = c(4, 3)
esti_mv = optim(inicial, ll, x = muestra, t = 0.5, hessian = TRUE)
sqrt(diag(solve(esti_mv$hessian)))

#### Creaci?n de la funci?n acumulada y adicionales ####

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
## Varianza ##
vtweibull <- function(qt, alpha, t) {
    if (t > 1 | t < 0) {stop("par?metro t solo puede contener valores entre 0 y 1")}
    if (qt < 0) {stop("par?metro qt solo puede ser superior a 0")}
    if (alpha < 0) {stop("par?metro alpha solo puede ser superior a 0")}
    ct = (-log(1 - t)) ^ (1 / alpha)
    var = ((qt ^ 2) / (ct) ^ (1 / alpha)) * (gamma(1 + (2 / alpha)) - (gamma(1 + (1 / alpha))) ^ 2)
    return(var)
}

## Valor esperado ##
mtweibull <- function(qt, alpha, t) {
    if (t > 1 | t < 0) {stop("par?metro t solo puede contener valores entre 0 y 1")}
    if (qt < 0) {stop("par?metro qt solo puede ser superior a 0")}
    if (alpha < 0) {stop("par?metro alpha solo puede ser superior a 0")}
    ct = (-log(1 - t)) ^ (1 / alpha)
    mt = ((qt) / (ct ^ (1 / alpha))) * gamma(1 + 1 / alpha)
    return(mt)
}

qtseq = seq(0, 10, 0.1)

#### Gráfico de valor esperado ####

e_x1 = data.frame(a = mtweibull(qtseq, 1, 0.5))
e_x1 = cbind(e_x1, data.frame(b = mtweibull(qtseq, 2, 0.5)))
e_x1 = cbind(e_x1, data.frame(c = mtweibull(qtseq, 3, 0.5)))
e_x1 = cbind(e_x1, data.frame(d = mtweibull(qtseq, 4, 0.5)))
e_x1 = cbind(e_x1, qtseq)

e_x2 = data.frame(a = mtweibull(1, qtseq, 0.5))
e_x2 = cbind(e_x2, data.frame(b = mtweibull(2, qtseq, 0.5)))
e_x2 = cbind(e_x2, data.frame(c = mtweibull(3, qtseq, 0.5)))
e_x2 = cbind(e_x2, data.frame(d = mtweibull(4, qtseq, 0.5)))
e_x2 = cbind(e_x2, qtseq)

## Valor esperado; eje horizontal es qt

esp1 <- ggplot(e_x1, aes(x=qtseq)) +
    geom_line(aes(y = a, colour = "#002b36")) + #alpha =1
    geom_line(aes(y = b, colour = "#268bd2")) + #alpha =2
    geom_line(aes(y = c, colour = "#859900")) + #alpha =3
    geom_line(aes(y = d, colour = "#dc322f")) + # alpha=4
scale_color_discrete(name="Valores de alpha",labels=c("1","2","3","4")) +
    labs(title="Valor esperado para Y~W(t=0.5)", 
         x="Parámetro q_t",y="Valor esperado",
         subtitle = "Valores fluctuantes del parámetro alpha") +
    ggthemes::theme_stata(scheme = "s1color",base_size = 10) + 
    theme(legend.position = "bottom")
## Valor esperado; eje horizontal es alpha

esp2 <- ggplot(e_x2, aes(x=qtseq)) +
    geom_line(aes(y = a, colour = "#002b36"), size = 1) + #qt =1
    geom_line(aes(y = b, colour = "#268bd2"), size = 1) + #qt =2
    geom_line(aes(y = c, colour = "#859900"), size = 1) + #qt =3
    geom_line(aes(y = d, colour = "#dc322f"), size = 1) + #qt= 4
    scale_y_continuous(limits=c(0,15)) + 
    scale_x_continuous(limits = c(0.50,10)) +#qt= 4
    scale_color_discrete(name="Valores de qt",labels=c("1","2","3","4")) +
    labs(title="Valor esperado para Y~W(t=0.5)", 
         x="Parámetro alpha",y="Valor esperado",
         subtitle = "Valores fluctuantes del parámetro qt") +
    ggthemes::theme_stata(scheme = "s1color",base_size = 10) + 
    theme(legend.position = "bottom")

plot3<- ggarrange(esp1,esp2,ncol=2,nrow=1)
plot3

#### Gráfico de varianza ####
## Varianza; eje horizontal es qt ##
v_x1 = data.frame(a = vtweibull(qtseq, 1, 0.5))
v_x1 = cbind(v_x1, data.frame(b = vtweibull(qtseq, 2, 0.5)))
v_x1 = cbind(v_x1, data.frame(c = vtweibull(qtseq, 3, 0.5)))
v_x1 = cbind(v_x1, data.frame(d = vtweibull(qtseq, 4, 0.5)))
v_x1 = cbind(v_x1, qtseq)

## Varianza; eje horizontal es alpha
v_x2 = data.frame(a = vtweibull(1, qtseq, 0.5))
v_x2 = cbind(v_x2, data.frame(b = vtweibull(2, qtseq, 0.5)))
v_x2 = cbind(v_x2, data.frame(c = vtweibull(3, qtseq, 0.5)))
v_x2 = cbind(v_x2, data.frame(d = vtweibull(4, qtseq, 0.5)))
v_x2 = cbind(v_x2, qtseq)

## Plot de varianza; eje horizontal es qt ##

var1 <- ggplot(v_x1,aes(x = qtseq)) +
    geom_line(aes(y = a, colour = "#002b36")) + #alpha =1
    geom_line(aes(y = b, colour = "#268bd2")) + #alpha =2
    geom_line(aes(y = c, colour = "#859900")) + #alpha =3
    geom_line(aes(y = d, colour = "#dc322f")) + # alpha = 4
    scale_color_discrete(name="Valores de alpha",labels=c("1","2","3","4")) +
    labs(title="Varianza para Y~W(t=0.5)", 
         x="Parámetro q_t",y="Varianza",
         subtitle = "Valores fluctuantes del parámetro alpha") +
    ggthemes::theme_stata(scheme = "s1color",base_size = 10) + 
    theme(legend.position = "bottom")

## Varianza; eje horizontal es alpha

var2 <- ggplot(v_x2) +
    geom_line(aes(x = qtseq, y = a, colour = "#002b36"), size = 1) + #qt =1
    geom_line(aes(x = qtseq, y = b, colour = "#268bd2"), size = 1) + #qt =2
    geom_line(aes(x = qtseq, y = c, colour = "#859900"), size = 1) + #qt =3
    geom_line(aes(x = qtseq, y = d, colour = "#dc322f"), size = 1) + 
    scale_y_continuous(limits=c(0,10)) + 
    scale_x_continuous(limits = c(1,10)) +#qt= 4
    scale_color_discrete(name="Valores de qt",labels=c("1","2","3","4")) +
    labs(title="Varianza para Y~W(t=0.5)", 
         x="Parámetro alpha",y="Varianza",
         subtitle = "Valores fluctuantes del parámetro qt") +
    ggthemes::theme_stata(scheme = "s1color",base_size = 10) + 
    theme(legend.position = "bottom")

plot2<- ggarrange(var1,var2,ncol=2,nrow=1)
plot2
#### Plot de densidad ####
qrv = data.frame(a = dtweibull(qtseq, qt = 1, alpha = 2, t = 0.5))
qrv = cbind(qrv, data.frame(b = dtweibull(qtseq, qt = 2, alpha = 2, t = 0.5)))
qrv = cbind(qrv, data.frame(c = dtweibull(qtseq, qt = 3, alpha = 2, t = 0.5)))
qrv = cbind(qrv, data.frame(d = dtweibull(qtseq, qt = 4, alpha = 2, t = 0.5)))
qrv = cbind(qrv, data.frame(e = dtweibull(qtseq, qt = 5, alpha = 2, t = 0.5)))

arv = data.frame(a = dtweibull(qtseq, qt = 3, alpha = 1, t = 0.5))
arv = cbind(arv, data.frame(b = dtweibull(qtseq, qt = 3, alpha = 2, t = 0.5)))
arv = cbind(arv, data.frame(c = dtweibull(qtseq, qt = 3, alpha = 3, t = 0.5)))
arv = cbind(arv, data.frame(d = dtweibull(qtseq, qt = 3, alpha = 4, t = 0.5)))
arv = cbind(arv, data.frame(e = dtweibull(qtseq, qt = 3, alpha = 5, t = 0.5)))

## Primer plot de densidad, qt=3, t=0.5

density1 <- ggplot(arv, aes(x=qtseq)) +
    geom_line(aes(y = a, colour = "#002b36")) + #alpha =1
    geom_line(aes(y = b, colour = "#268bd2")) + #alpha = 2
    geom_line(aes(y = c, colour = "#859900")) + #alpha = 3
    geom_line(aes(y = d, colour = "#dc322f")) + # alpha=4
    geom_line(aes(y = e, colour = "#b58900")) +
    scale_color_discrete(name="Valores de alpha",labels=c("1","2","3","4","5")) +
    labs(title="Función de densidad para Y~W(q_{t}=3,t=0.5)", 
         x="Y",y="Probabilidad",
         subtitle = "Valores fluctuantes del parámetro de alpha") + 
    scale_x_continuous(limits = c(0,10))+scale_y_continuous(limits=c(0,1)) +
    ggthemes::theme_stata(scheme = "s1color",base_size = 10) + 
    theme(legend.position = "bottom")

## Segundo plot de densidad, alpha = 5, t=0.5
density2 <- ggplot(qrv) +
    geom_line(aes(y = a, x = qtseq, colour = "#002b36")) + #qt =1
    geom_line(aes(y = b, x = qtseq, colour = "#268bd2")) + #qt = 2
    geom_line(aes(y = c, x = qtseq, colour = "#859900")) + #qt = 3
    geom_line(aes(y = d, x = qtseq, colour = "#dc322f")) + #qt=4
    geom_line(aes(y = e, x = qtseq, colour = "#b58900")) + #qt=5    
    scale_color_discrete(name="Valores de qt",labels=c("1","2","3","4","5")) +
    labs(title="Función de densidad para Y~W(alpha = 2, t=0.5)", x="Y",y="Probabilidad",subtitle= "Valores fluctuantes del parámetro qt" ) + 
    scale_y_continuous(limits=c(0,1)) + 
    ggthemes::theme_stata(scheme = "s1color",base_size = 10) + 
    theme(legend.position = "bottom")

plot1 <- ggarrange(density1,density2,ncol=2,nrow=1)
plot1

#### Segunda parte: creaci?n de la funci?n de regresi?n #####

logptweibull <- function(value, qt, sigma, t) {
    # if (t > 1 | t < 0) { stop("par?metro t solo puede contener valores entre 0 y 1") }
    # if (qt <= 0) { stop("par?metro qt solo puede ser superior a 0") }
    # if (sigma <= 0) { stop("par?metro sigma solo puede ser superior a 0") }
    alpha = 1 / (log(sigma + 1))
    ct = (-log(1 - t)) ^ (1 / alpha)
    b = qt / ct
    a = alpha
    pt = pweibull(value, a, b,log=TRUE)
    return(pt)
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
info <- lancet()

reg_ces_wei = function(data, li, lf, t) {
    ### Test ###
    # data = VADeaths
    # li = 2
    # li = 3
    # param = c(rep(0.5, ncol(data) + 1))
    # t=0.5
    data_l_inf = data[,li]
    data_l_sup = data[, lf]
    covar = as.matrix(data[, - c(li, lf)])
    n = ncol(covar)
    d = nrow(covar)

    ll = function(param, l_inf,l_fin, t) {
        ### Creaci?n de betas ###
        b_0 = 0
        b = as.matrix(rep(0, ncol(covar)))
        ### Asignaci?n de betas ###
        b_0 = param[1]        
        for (i in 1:n) {b[i] = param[i + 1]}
        ### Asignaci?n de sigma ###
        sigma = param[length(param)]
        ### Asignaci?n de ct ###
        ct = -(log(1-t))

        ### Creaci?n de funci?n de enlace ###
        qt = exp(b_0 + covar %*% b)
        ### Inicializaci?n de la log-verosimilitud ###
        s_ll = 0

        ### C?lculo de la log-verosimilitud s###
        for (j in 1:d) {
        qt_j = qt[j,]
        
        #s_ll = s_ll - logDiffExp(logptweibull(data_l_sup[j], qt_j, sigma, t),logptweibull(data_l_inf[j], qt_j, sigma, t))
        s_ll = s_ll - log(ptweibull(data_l_sup[j], qt_j, sigma, t)- ptweibull(data_l_inf[j], qt_j, sigma, t))
        }
        return(s_ll)
    }
    
    #t_index = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
    #inicial= list(c(7.6601709,-0.2040200, 0.3173124),
               #c(7.866936,-0.204009,0.317338),
               #c(7.9962011,-0.2040046,0.3173342),
               #c(8.0952017,-0.2040035, 0.3173330),
               #c(8.1793198,-0.2040044,0.3173373),
               #c(8.2562373,-0.2040037,0.3173392),
               #c(8.3314962,-0.2040040, 0.3173363),
               #c(8.4114690,-0.2039752, 0.3173351),
               #c(8.5102048,-0.2040047, 0.3173366))
    
    inicial = c(as.vector(coef(lm((data_l_inf + data_l_sup)/2~covar,data=data.frame(cbind(data_l_inf,covar))))),1)
    #fit_mv =   optim(par = inicial, method="L-BFGS-B",fn = ll,hessian = T, t = t, upper = c(rep(Inf, ncol(covar) + 2)),lower = c(rep(0, ncol(covar) + 1), 0.1))
    fit_mv = nlminb(inicial,ll,l_inf = li, l_fin = lf,t = t,lower = c(rep(0,ncol(covar)+1,0.1)),upper =c(rep(Inf,ncol(covar)+2)))
    inicial = fit_mv$par
    return(fit_mv)
}


#te = seq(0.1,0.9,0.1)
#fit_mv_0.1 = reg_ces_wei(info, 2, 3, 0.1)
#fit_mv_0.2 = reg_ces_wei(info, 2, 3, 0.2)
#fit_mv_0.3 = reg_ces_wei(info, 2, 3, 0.3)
#fit_mv_0.4 = reg_ces_wei(info, 2, 3, 0.4)
#fit_mv_0.5 = reg_ces_wei(info, 2, 3, 0.5)
#fit_mv_0.6 = reg_ces_wei(info, 2, 3, 0.6)
#fit_mv_0.7 = reg_ces_wei(info, 2, 3, 0.7)
#fit_mv_0.8 = reg_ces_wei(info, 2, 3, 0.8)
#fit_mv_0.9 = reg_ces_wei(info, 2, 3, 0.9)

#for (i in 1:length(te)) {print(eval(parse(text = paste("fit_mv_", te[i], "$par", sep = ""))))}
#conf.level <- 0.95
#crit <- qnorm((1+conf.level)/2)





head(sim)


fit_mv_0.5$par[1] + c(-1,1)*crit*sqrt(diag(solve(fit_mv_0.5$hessian))[1])
fit_mv_0.5$par[2] + c(-1,1)*crit*sqrt(diag(solve(fit_mv_0.5$hessian))[2])
fit_mv_0.5$par[3] + c(-1,1)*crit*sqrt(diag(solve(fit_mv_0.5$hessian))[3])
fit_mv_0.5$par[4] + c(-1,1)*crit*sqrt(diag(solve(fit_mv_0.5$hessian))[4])
fit_mv_0.5$par[5] + c(-1,1)*crit*sqrt(diag(solve(fit_mv_0.5$hessian))[5])

sesgo = sum(fit_mv_0.5$par - c(0,b_1,b_2,b_3,sig));sesgo
ecm = sum(fit_mv_0.5$par - c(0,b_1,b_2,b_3,sig))**2;ecm

M = 100
t_sim = 0.5

fit_sim_0.1= c()
fit_sim_0.2= c()
fit_sim_0.3= c()
fit_sim_0.4= c()
fit_sim_0.5= c()
fit_sim_0.6= c()
fit_sim_0.7= c()
fit_sim_0.8= c()
fit_sim_0.9= c()

for(j in 1:length(t_sim)){
    X2 = rbeta(10000,2,3)
    X3 = rnorm(10000,3,0.5)
    X4 = rgamma(10000,2,4)

    b_1 = 0.3
    b_2 = 0.6
    b_3 = 0.8
    sig = 2
    Qt  = exp(b_1*X2 + b_2*X3 + b_3*X4)
    Yh = c()

fit_sim_0.5 = list()
   for(k in 1:M){
        for(i in 1:length(Qt)){
            Yh[i]= rtweibull(1,Qt[i],sig,t_sim[j])
            }
        sim = data.frame(   Li = ifelse(Yh < quantile(Yh,0.25),0,
                                ifelse(Yh <  quantile(Yh,0.5),quantile(Yh,0.25),
                                ifelse(Yh < quantile(Yh,0.75),quantile(Yh,0.50),quantile(Yh,0.70)))),
                                Ls = ifelse(Yh < quantile(Yh,0.25),quantile(Yh,0.25),
                                ifelse(Yh <  quantile(Yh,0.5),quantile(Yh,0.5),
                                ifelse(Yh < quantile(Yh,0.75),quantile(Yh,0.75),quantile(Yh,0.99)))),
                                X2 =X2,
                                X3 = X3,
                                X4 = X4)
        print(k)
        fit_sim_0.5  <- append(fit_sim_0.5, reg_ces_wei(sim,1,2,t_sim[j]))
    }
}

fit_sim_0.5