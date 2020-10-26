real_data <- lancet()

t_sim = seq(0.25,0.9,0.25)
var = list()
for (j in 1:length(t_sim)) {
  var = append(var,list(reg_ces_wei(data = real_data,5,6,t_sim[j])))
}

val <- matrix(ncol=length(var[[1]]$par),nrow = 0)
for (l in 1:length(t_sim)) {
    pba <- var[[l]]
    val <- rbind(val,pba$par)
}
val

for (i in 1:length(t_sim)){
  print(round(coef(rq(li~sexo+C2P9+C2P27+C2P21,real_data,tau = t_sim[i])),5))
}

