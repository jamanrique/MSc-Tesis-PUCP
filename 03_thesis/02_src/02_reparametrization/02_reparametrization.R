
#### Creaci?n de la funci?n de distribuci?n para n?meros aleatorios ####

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
