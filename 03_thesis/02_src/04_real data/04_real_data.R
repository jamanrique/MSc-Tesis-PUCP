real_data <- lancet()

t_sim = seq(0.1,0.9,0.1)
var = list()
for (j in 1:length(t_sim)) {
  var = append(var,list(reg_ces_wei(data = real_data,5,6,t_sim[j])))
}
var
val <- matrix(ncol=length(var[[1]]$par),nrow = 0)
for (l in 1:length(t_sim)) {
    pba <- var[[l]]
    val <- rbind(val,pba$par)
}
val

library(gamlss)
library(gamlss.cens)

gen.cens(family="WEI3",type="interval")
gamlss(Surv(li,lf,type="interval2")~.,family = WEI3ic,data = real_data)

fisher_info <- solve(var[[5]]$hessian)
prop_sigma <- sqrt(diag(fisher_info)); prop_sigma <- diag(prop_sigma);prop_sigma
limit_upper <- var[[5]]$par + 1.96*prop_sigma; limit_upper
limit_lower <- var[[5]]$par - 1.96*prop_sigma; limit_lower

rbind(limit_upper,limit_lower)
