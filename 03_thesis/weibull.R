library(rgl)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(latex2exp)
library(gridExtra)

rm(list=ls())
x <- 10000
seqsigma <- seq(1,9,1)
col_names <- c("\U003C3 = 1","\U003C3 = 2","\U003C3 = 3","\U003C3 = 4","\U003C3 = 5","\U003C3 = 6","\U003C3 = 7","\U003C3 = 8","\U003C3 = 9")

dbf <- data.frame(rweibull(x,2,1),"\U003C3 = 1")
for (i in 2:length(seqsigma)) {
  dbf_aux <- data.frame(rweibull(x,2,seqsigma[i]),col_names[i])
  dbf <- rbind(dbf, setNames(dbf_aux,names(dbf)))
}
colnames(dbf) <- c("value","x1")

density1 <- ggplot(dbf, aes(x=value,y=x1,fill=x1)) +
  geom_density_ridges(aes(point_fill=x1),  alpha = .35, point_alpha = 1) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_fill", values = c(21, 22, 23,24,25,26,27,28,29)) + 
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(axis.text.y = element_blank()) +
  labs(x="y",y="Densidad",title="Y \U0007E W(\U03B1,\U003C3), \U03B1 = 2",subtitle = "Simulación de 10,000 valores") +
  theme(plot.title = element_text(face="bold",hjust=1,size=16),
        plot.subtitle = element_text(face="italic",hjust=1,size=14),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "gray40", size = .25),
        panel.grid.minor = element_line(color = "gray70", size = .25),
        axis.title = element_text(size=8))


seqalpha <- seq(1,9,1)

col_names <- c("\U03B1 = 1","\U03B1 = 2","\U03B1 = 3","\U03B1 = 4","\U03B1 = 5","\U03B1 = 6","\U03B1 = 7","\U03B1 = 8","\U03B1 = 9")

dbf_b <- data.frame(data.frame(rweibull(x,1,2)),"\U03B1 = 1")
for (i in 2:length(seqalpha)) {
  dbf_aux <- data.frame(rweibull(x,seqalpha[i],2),col_names[i])
  dbf_b <- rbind(dbf_b, setNames(dbf_aux,names(dbf_b)))
}

colnames(dbf_b) <- c("value","x1")

density2 <- ggplot(dbf_b, aes(x=value,y=x1,fill=x1)) +
  geom_density_ridges(aes(point_fill=x1),  alpha = .35, point_alpha = 1) +
  scale_point_color_hue(l = 40) +
  scale_discrete_manual(aesthetics = "point_fill", values = c(21, 22, 23,24,25,26,27,28,29)) +
  theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  theme(axis.text.y = element_blank()) +
  labs(x="y",y="Densidad" , title="Y \U0007E W(\U03B1,\U003C3), \U003C3 = 2",subtitle = "Simulación de 10,000 valores")+
  theme(plot.title = element_text(face="bold",hjust=1,size=16),
        plot.subtitle = element_text(face="italic",hjust=1,size=14),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.major = element_line(color = "gray40", size = .25),
        panel.grid.minor = element_line(color = "gray70", size = .25),
        axis.title = element_text(size=8))
 

ggarrange(density1,density2,ncol=2,legend="bottom")

#### Generación de densidad ####

rm(list=ls())

rtweibull = function(n, qt, alpha, t) {
  ## Reparametrizaci?n de los datos a la f?rmula b?sica de Weibull en R ##
  beta = qt / (-log(1 - t)) ^ (1 / alpha)
  qm = rweibull(n = n, scale = beta, shape = alpha)
  return(qm)
}

seq_qt <- seq(1,10,length.out = 5)
seq_tau <- seq(0.1,0.9,length.out = 5)
alpha <- 1
n <- 10000

densityplot <- function(n,qt,tau,alpha,xlim1,xlim2){
  dens <- data.frame(x=rtweibull(n = n,qt = qt,alpha = alpha,t = tau))
  dens.mean <- data.frame(y=mean(dens$x))
  gg <- ggplot(dens,
               aes(x=x)) +
    geom_density(color="black",fill="firebrick",alpha=0.8) +
    xlim(c(xlim1,xlim2))+
    labs(y="Densidad",
         x = "y", 
         title = (paste("Y \U0007E W\U01D63(Qt = ",qt,", \U003B1 = ",alpha,"); \U003C4 =",tau,sep = "")),
         subtitle = "Simulación de 10,000 valores")+
    theme(plot.title = element_text(face="bold",hjust=1,size=10),
          plot.subtitle = element_text(face="italic",hjust=1,size=8),
          panel.background = element_rect(fill = "gray90"),
          panel.grid.major = element_line(color = "gray40", size = .25),
          panel.grid.minor = element_line(color = "gray70", size = .25),
          axis.title = element_text(size=8))
  return(gg)
}
gg_plots <- list()
for (qt in 1:length(seq_qt)) {
    for(t in 1:length(seq_tau)){
    gg_plots <- append(gg_plots,list(densityplot(n,seq_qt[qt],seq_tau[t],2,0,30)))
  }
}
plot_grid <- do.call(grid.arrange,gg_plots)

seq_alpha <- seq(1,10,length.out = 5)
gg_plots_a <- list()
for (qt in 1:length(seq_alpha)) {
  for(t in 1:length(seq_tau)){
    gg_plots_a <- append(gg_plots_a,list(densityplot(n,2,seq_tau[t],seq_alpha[qt],0,20)))
  }
}
plot_grid_alpha <- do.call(grid.arrange,gg_plots_a)

vtweibull <- function(qt, alpha, t) {
  ct = (-log(1 - t)) ^ (1 / alpha)
  var = ((qt^2) / (ct)) * (gamma(1 + (2 / alpha)) - (gamma(1 + (1 / alpha))) ^ 2)
  return(var)
}

mtweibull <- function(qt, alpha, t) {
  ct = (-log(1-t))^(1/alpha)
  mt = ((qt) / (ct)) * gamma(1 + 1 / alpha)
  return(mt)
}

mtweibull(10,2,0.5)
vtweibull(10,2,0.5)
