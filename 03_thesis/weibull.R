library(rgl)
library(tidyverse)
library(ggridges)

rm(list=ls())
x <- 1000
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
  labs(x="y",y="Función de densidad para Y \U0007E W(\U03B1,\U003C3), \U03B1 = 2") 


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
  labs(x="y",y="Función de densidad para Y \U0007E W(\U03B1,\U003C3), \U003C3 = 2") 

density1
density2
