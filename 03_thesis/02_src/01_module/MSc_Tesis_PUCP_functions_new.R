library(gamlss)
library(gamlss.cens)
gen.cens(family="WEI3",type="interval")

Ct = function(alpha, tau){
  (-log(1-tau))^(1/alpha)
}
F_Wr = function(Y,Qt,alpha,tau){
  B = Qt/Ct(alpha,tau)
  pweibull(Y,alpha,B)
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
  rweibull(n,alpha,B)
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
  min_Y = min(Y); Q8_Y = quantile(Y,0.8); interval = (Q8_Y-min_Y) / 6
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
  fit_mv = nlminb(start=param,objective=likelihood,betas=Bs,df=df,tau=tau)
  return(fit_mv)
}

sim = Simulation(10000)
betas_sim = c(7,0.3,0.84,2.5)
tau_sim = 0.2
alpha_sim = 2
df_sim = DF_Simulation(sim,betas_sim,alpha_sim,tau_sim)
df_sim = na.omit(df_sim)

m0 = gamlss(Surv(Li,Ls,type="interval2")~.,family = WEI3ic,data = na.omit(df_sim))
init = as.vector(c(coef(m0),m0$sigma.coefficients));init

reg_Wr(data = df_sim,li = 1,lf = 2,tau = tau_sim,param = init)
