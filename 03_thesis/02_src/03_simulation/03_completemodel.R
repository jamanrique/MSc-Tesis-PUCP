    # Justo Manrique - 20091107
    # Regresi?n cuant?lica en datos intervalares: un estudio sobre la disparidad de ingresos entre hombres y mujeres en el sistema de salud peruano
    
    #### Carga de librer?as ####
    
    #### Segunda parte: creaci?n de la funci?n de regresi?n #####
    
    fit_sim= list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    
    
    t_sim = seq(0.10,0.9,0.1)
    b_0 = 0.5
    b_1 = 0.3
    b_2 = 0.6
    b_3 = 0.8
    sig = 2
    
    ## incluir intercepto
    X2 = rbeta(10000,2,3)
    X3 = rnorm(10000,2,0.25)
    X4 = rgamma(10000,2,25)
    Qt  = exp(b_0+b_1*X2 + b_2*X3 + b_3*X4) ### vector 10000
    
    for(j in 1:length(t_sim)){
        Yh = c()
        M = 100
        while (M != 0) {
            try({
                for(i in 1:length(Qt)){Yh[i]= rtweibull(1,Qt[i],sig,t_sim[j])}
                intervals=7
                seq_interv= c(round(seq(min(Yh),quantile(Yh,0.80),(quantile(Yh,0.80)-min(Yh))/intervals),2),Inf)
                Ls = c()
                for (u in 1:length(Yh)) {
                    for (n in 1:length(seq_interv)) {
                        if (Yh[u] <seq_interv[n]) {
                            Ls[u] = seq_interv[n]
                            break
                        }}}
                Li = c()
                for (p in 1:length(Yh)) {
                    for (w in 1:length(seq_interv)) {
                        if (Yh[p] > rev(seq_interv)[w]) {
                            Li[p] = rev(seq_interv)[w]
                            break
                        }}}
                sim = data.frame(
                    Li = Li,
                    Ls = Ls,
                    X2 =X2,
                    X3 = X3,
                    X4 = X4)
                fit_sim[[j]]  <- append(fit_sim[[j]], list(reg_ces_wei(sim,1,2,t_sim[j])))
                print(paste("Resultado final de la simulación Nª ",abs(100-M+1),", para el cuantil ",t_sim[j], ": CONVIRGIÓ", sep = ""))
                M <- M - 1
                }
            )
        }
        M = 100
    }
    
    ## Hacer una regresión con del dato directo
    ## intervalos igualamente espaciados
    
    crit=qnorm(0.975)
    
    coverage = list(
        t_0.10 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.20 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.30 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.40 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.50 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.60 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.70 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.80 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list()),
        t_0.90 = list( b_0 = list(),b_1 = list(),b_2=list(),b_3=list(),sig=list())
        )
    
    sesgo_b1 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    sesgo_b2 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    sesgo_b3 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    sesgo_b0 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    sesgo_sig = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    
    ecm_b0 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    ecm_b1 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    ecm_b2 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    ecm_b3 = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    ecm_sig = list(
        t_0.10 = list(),
        t_0.20 = list(),
        t_0.30 = list(),
        t_0.40 = list(),
        t_0.50 = list(),
        t_0.60 = list(),
        t_0.70 = list(),
        t_0.80 = list(),
        t_0.90 = list()
    )
    
    for(j in 1:length(t_sim)){
        for(l in 1:M){
            pba = fit_sim[[j]][[l]]
            b_0_test = pba$par[1] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[1])
            b_1_test = pba$par[2] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[2])
            b_2_test = pba$par[3] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[3])
            b_3_test = pba$par[4] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[4])
            sig_test = pba$par[5] + c(-1,1)*crit*sqrt(diag(solve(pba$hessian))[5])
        
            if (between(b_0,b_0_test[1],b_0_test[2])) {coverage[[j]][[1]] = append(coverage[[j]][[1]],1)}
            if (between(b_1,b_1_test[1],b_1_test[2])) {coverage[[j]][[2]] = append(coverage[[j]][[2]],1)}
            if (between(b_2,b_2_test[1],b_2_test[2])) {coverage[[j]][[3]] = append(coverage[[j]][[3]],1)}
            if (between(b_3,b_3_test[1],b_3_test[2])) {coverage[[j]][[4]] = append(coverage[[j]][[4]],1)}
            if (between(sig,sig_test[1],sig_test[2])) {coverage[[j]][[5]] = append(coverage[[j]][[5]],1)}
            
            sesgo_b0[[j]] = append(sesgo_b0[[j]],sum(pba$par[1] - c(b_0)))
            sesgo_b1[[j]] = append(sesgo_b1[[j]],sum(pba$par[2] - c(b_1)))
            sesgo_b2[[j]] = append(sesgo_b2[[j]],sum(pba$par[3] - c(b_2)))
            sesgo_b3[[j]] = append(sesgo_b3[[j]],sum(pba$par[4] - c(b_3)))
            sesgo_sig[[j]] = append(sesgo_b3[[j]],sum(pba$par[5] - c(sig)))
    
            ecm_b0[[j]] = append(ecm_b0[[j]],sum(pba$par[1] - c(b_0))**2)
            ecm_b1[[j]] = append(ecm_b1[[j]],sum(pba$par[2] - c(b_1))**2)
            ecm_b2[[j]] = append(ecm_b1[[j]],sum(pba$par[3] - c(b_2))**2)
            ecm_b3[[j]] = append(ecm_b1[[j]],sum(pba$par[4] - c(b_3))**2)
            ecm_sig[[j]] = append(ecm_b1[[j]],sum(pba$par[5] - c(sig))**2)
    }
    }
    
t_cuadro <- data.frame()

for(p in 1:length(t_sim)){
    sr_b0 <- mean(unlist(sesgo_b0[[p]]))/b_0
    sr_b1 <- mean(unlist(sesgo_b1[[p]]))/b_1
    sr_b2 <- mean(unlist(sesgo_b2[[p]]))/b_2
    sr_b3 <- mean(unlist(sesgo_b3[[p]]))/b_3
    sr_sig <- mean(unlist(sesgo_sig[[p]]))/sig
    cov_b0 <-sum(unlist(coverage[[p]][[1]]))/100
    cov_b1 <-sum(unlist(coverage[[p]][[2]]))/100
    cov_b2 <-sum(unlist(coverage[[p]][[3]]))/100
    cov_b3 <-sum(unlist(coverage[[p]][[4]]))/100
    cov_sig <-sum(unlist(coverage[[p]][[5]]))/100
    ecmb0 <- mean(unlist(ecm_b0[[p]]))
    ecmb1 <- mean(unlist(ecm_b1[[p]]))
    ecmb2 <- mean(unlist(ecm_b2[[p]]))
    ecmb3 <- mean(unlist(ecm_b3[[p]]))
    ecmsig <- mean(unlist(ecm_sig[[p]]))
    
    tot_sr <- rbind(sr_b0,sr_b1,sr_b2,sr_b3,sr_sig)
    tot_cov <- rbind(cov_b0,cov_b1,cov_b2,cov_b3,cov_sig)
    tot_ecm <- round(rbind(ecmb0,ecmb1,ecmb2,ecmb3,ecmsig),4)
    cuad <- as.data.frame(cbind(tot_sr,tot_cov,tot_ecm))
    t_cuadro <- rbind(t_cuadro,cuad);t_cuadro
}
    
val <- c()
for (i in 1:100) {
    pba <- fit_sim$t_0.10[[i]]
    val <- append(val,pba$value)
    
}
plot(density(val))
    
