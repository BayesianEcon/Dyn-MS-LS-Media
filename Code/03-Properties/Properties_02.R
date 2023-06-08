
list.of.packages <- c("ggplot2", "ggpubr", "patchwork", "tidyr", "scales", "igraph")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")

library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyr)
library(scales)
library(igraph)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

start_time <- Sys.time()  

avgDeg = function(n, alpha, d, sigma2){
  
  res = (n-1)*exp(alpha)*(4*sigma2+1)^(-d/2)
  return(res)
}

N = 100
sigma = seq(0.01,3.01, by = 0.01)
alpha = seq(0.01,3.01, by = 0.01)
mean_edge<-rep(0,length(sigma))
zeta2<- rep(0,N)

Alist<-list()

for(i in 1:length(sigma)){
  zeta1<- rnorm(N, 0, sqrt(sigma[i]))
  #zeta2<- rnorm(N, 0, sigma[i])
  D<-as.matrix(dist(cbind(zeta1, zeta2)))
  A = exp(2 - D^2)
  A = apply(A, MARGIN= c(1,2), FUN = function(x){rpois(1,x)})
  diag(A)<-0
  
  mean_edge[i]<-  mean(apply(A, 1, FUN = sum ))
  print(i)
}


for(i in 1:length(sigma)){
  zeta1<- rnorm(N, 0, sqrt(sigma[i]))
  #zeta2<- rnorm(N, 0, sigma[i])
  D<-as.matrix(dist(cbind(zeta1, zeta2)))
  A = exp(-0.1 - D^2)
  A = apply(A, MARGIN= c(1,2), FUN = function(x){median(rbinom(10,1,x))})
  diag(A)<-0
  
  mean_edge[i]<-  mean(apply(A, 1, FUN = sum ))
  print(i)
}


for(i in 1:length(sigma)){
  zeta1<- rnorm(N, 0, 1)
  zeta2<- rep(0,N)
  #zeta2<- rnorm(N, 0, 1)
  D<-as.matrix(dist(cbind(zeta1, zeta2)))
  A = exp(alpha[i] - D^2)
  A= A[upper.tri(A)]#A = exp(alpha[i] - D^2)
  A = sapply(A, FUN = function(x){rpois(1,x)})


  mean_edge[i]<-mean(A)
  print(i)
}


N = 100
sigma2 = seq(0.01,3.01, by = 0.02)
alpha = seq(0.01, 3.01, by = 0.02)

mean_edge_m<-matrix(0,length(alpha), length(sigma2)+1)
sd_edge_m<-matrix(0, length(alpha),length(sigma2)+1)
trans<-matrix(0, length(alpha),length(sigma2)+1)

zeta2<- rep(0, N)


for(i in 1:length(alpha)){
  

  alf<-alpha[i]
  
  for(j in 0:length(sigma2)){
    
  if(j ==0){
    
    zeta1<- rep(0, N)
   # zeta1<- rnorm(N, 0, sqrt(sigma2[j]))
    
    D<-as.matrix(dist(cbind(zeta1, zeta2)))
    A = exp(alf - D^2)
    A = apply(A, MARGIN= c(1,2), FUN = function(x){(rpois(1,x))})
    diag(A)<-0
    
    mean_edge_m[i,j+1]<-  mean(apply(A, 1, FUN = sum ))
    sd_edge_m[i,j+1]<-  sd(apply(A, 1, FUN = sum ))
    
    g= graph_from_adjacency_matrix(A, mode = "undirected", weighted = T)
    tr<-mean(transitivity(g, type = "barrat"), na.rm= T)
    
    trans[i, j+1]<-tr
    
  }else{
    
    
    zeta1<- rnorm(N, 0, sqrt(sigma2[j]))
    #zeta2<- rnorm(N, 0, sqrt(sigma2[j]))
    
    D<-as.matrix(dist(cbind(zeta1, zeta2)))
    A = exp(alf - D^2)
    A = apply(A, MARGIN= c(1,2), FUN = function(x){(rpois(1,x))})
    diag(A)<-0
    
    g= graph_from_adjacency_matrix(A, mode = "undirected", weighted = T)
    tr<-mean(transitivity(g, type = "barrat"), na.rm= T)
    
    
    mean_edge_m[i,j+1]<-  mean(apply(A, 1, FUN = sum ))
    sd_edge_m[i,j+1]<-  sd(apply(A, 1, FUN = sum ))
    
    trans[i, j+1]<-tr
    
  }
    
  }
  
  print(i)
}

colnames(mean_edge_m)<- c(0,sigma2)
rownames(mean_edge_m)<- alpha

colnames(sd_edge_m)<- c(0,sigma2)
rownames(sd_edge_m)<- alpha

colnames(trans)<- c(0,sigma2)
rownames(trans)<- alpha

DB<-reshape2::melt(mean_edge_m)
library(ggplot2)

colnames(DB)<-c("alpha", "sigma2", "value")

DB_sd<-reshape2::melt(sd_edge_m)
colnames(DB_sd)<-c("alpha", "sigma2", "value")


library(ggplot2)

DB_sd$di <-DB_sd$value^2/(DB$value)


avgDeg = function(n, alpha, d, sigma2){
  
  res =(n-1)*exp(alpha)*(4*sigma2+1)^(-d/2)
  return(res)
}



Var = function(n, alpha, d, sigma2){
  
  
  a = exp(2*alpha)*(n-1)*(8*sigma2+1)^(-d/2)+(n-1)*(n-2)*exp(2*alpha)*(2*sigma2+1)^(-d/2)*(6*sigma2+1)^(-d/2)
  b = (n-1)*exp(alpha)*(4*sigma2+1)^(-d/2)
  
  di =  a + b - b^2
  
  return(di)
  
}


DI = function(n, alpha, d, sigma2){
  
  # a = 2*(n-1)*exp(2*alpha)*(2*sigma2 + 0.5 )^(-d/2)*( (4/(4+sigma2)) + (1/sigma2) )^(-d/2)*(sigma2)^(-d/2)*2^(-d/2) +  (n-1)*(n-2)*exp(2*alpha)*(sigma2+0.5)^(-d)*( (2/(2+sigma2)) + (1/sigma2))^(-d/2)*(sigma2)^(-d/2)*2^(-d)
  a = exp(2*alpha)*(n-1)*(8*sigma2+1)^(-d/2)+(n-1)*(n-2)*exp(2*alpha)*(2*sigma2+1)^(-d/2)*(6*sigma2+1)^(-d/2)
  b = (n-1)*exp(alpha)*(4*sigma2+1)^(-d/2)
  
  di =  1+ (a/b) - b
  
  return(di)
  
}


n = 100

sigma2 = seq(0.01,3.01, by = 0.02)
alpha = seq(0.01, 3.01, by = 0.02)
             
Dbtrue<- data.frame( sigma2 = rep(sigma2, length(alpha)), alpha = rep(sigma2, each=length(sigma2)) )
Dbtrue$avg <- avgDeg(n=100, Dbtrue$alpha, 1, Dbtrue$sigma2)
Dbtrue$std <- sqrt(Var(n=100, Dbtrue$alpha, 1, Dbtrue$sigma2))
Dbtrue$dix <- DI(n=100, Dbtrue$alpha, 1, Dbtrue$sigma2)

end_time <- Sys.time()
end_time-start_time

alpha = seq(0.01, 2.01, 0.01)
sigma2 = seq(0.01, 3.01, 0.01)

p1<-ggplot(DB[DB$alpha %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= sigma2, y = value  , group = factor(alpha) ,color = alpha)) + labs(title = "Avg Strength - Sigma2", x =   expression(sigma^2), y = "Average Strength",color=expression(alpha), linetype=expression(alpha))+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p1<- p1 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")
p1<- p1 + geom_line(data=Dbtrue[Dbtrue$alpha %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x = sigma2, y = avg, linetype = factor(alpha)))

p2<-ggplot(DB[DB$sigma2 %in% c(0,0.01,0.11,0.51,1.01,1.1,1.51,2),], aes(x= alpha, y = value  ,  group = factor(sigma2), color = sigma2))+ labs(title = "Avg Strength - Alpha", x = expression(alpha) , y = "Average Strength" ,color=expression(sigma^2), linetype =expression(sigma^2))+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p2<- p2 + geom_line(data=Dbtrue[Dbtrue$sigma2 %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x = alpha, y = avg, linetype = factor(sigma2)))
p2<- p2 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")
library(patchwork)
p<-p1|p2
p

colnames(DB_sd)<-c("alpha", "sigma2", "value")

p3<-ggplot(DB_sd[DB_sd$alpha %in% c(0,0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= sigma2, y = value  , group = factor(alpha), color = alpha)) + labs(title = "SD Strength - Sigma2", x =   expression(sigma^2) , y = "SD Strength",color= expression(alpha), linetype= expression(alpha))+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p3<- p3 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")
p3<- p3 + geom_line(data=Dbtrue[Dbtrue$alpha %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x = sigma2, y = std, linetype = factor(alpha)))

p4<-ggplot(DB_sd[DB_sd$sigma2 %in% c(0,0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= alpha, y = value  , group = factor(sigma2),color = sigma2))+ labs(title = "SD Strength - Alpha",  x = expression(alpha),  y = "SD Strength",color= expression(sigma^2), linetype= expression(sigma^2) )+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p4<- p4 + geom_line(data=Dbtrue[Dbtrue$sigma2 %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x = alpha, y = std, linetype = factor(sigma2)))
p4<- p4 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")

library(patchwork)
p<-p3|p4
p

colnames(DB_sd)<-c("alpha", "sigma2", "value", "di")


p5<-ggplot(DB_sd[DB_sd$alpha %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= sigma2, y = di , group = factor(alpha) ,color = alpha)) + labs(title = "Dispersion Index - Strength - Sigma2" ,  x =   expression(sigma^2), y = "DI",color= expression(alpha), linetype=expression(alpha),)+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p5<-p5+geom_hline(yintercept=1, linetype = "dashed")
p5<- p5 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")
p5<- p5 + geom_line(data=Dbtrue[Dbtrue$alpha %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x = sigma2, y = dix, linetype = factor(alpha)))

p6<-ggplot(DB_sd[DB_sd$sigma2 %in% c(0,0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= alpha, y = di , group = factor(sigma2) ,color = sigma2))+ labs(title = "Dispersion Index - Strength - Alpha",  x = expression(alpha),y = "DI",color=expression(sigma^2), linetype=expression(sigma^2))+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p6<-p6+geom_hline(yintercept=1, linetype = "dashed")
p6<- p6 + geom_line(data=Dbtrue[Dbtrue$sigma2 %in% c(0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x = alpha, y = dix, linetype = factor(sigma2)))
p6<- p6 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")

library(patchwork)
p<-(p1|p2)/(p3|p4)/(p5|p6) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
p



#trans<-trans[,-1]
DB_trans<-reshape2::melt(trans)
colnames(DB_trans)<-c("alpha", "sigma2", "value")


p7<-ggplot(DB_trans[DB_trans$alpha %in% c(0,0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= sigma2, y = value , group = factor(alpha) ,color = alpha)) + labs(title = "Weighted CC (Barrat) - Sigma2" ,  x =   expression(sigma^2), y = "CC",color= expression(alpha), linetype= expression(alpha))+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p7<-p7+geom_hline(yintercept=1, linetype = "dashed")
p7<- p7 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")
p8<-ggplot(DB_trans[DB_trans$sigma2 %in% c(0,0.01,0.11,0.51,1.01,1.11,1.51,2),], aes(x= alpha, y = value , group = factor(sigma2) ,color = sigma2))+ labs(title = "Weighted CC (Barrat) - Alpha",  x = expression(alpha),y = "CC",color=expression(sigma^2), linetype= expression(sigma^2))+ geom_line() + theme_classic() + theme(legend.position = "bottom")
p8<-p8+geom_hline(yintercept=1, linetype = "dashed")
p8<- p8 + scale_color_gradient(low = "#d4d4d4",high = "#da6d50")


pp1<-(p1|p2)/(p3|p4) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
pp1

ggsave(pp1, filename = "Figures/Properties/FigureB.1.pdf", width = 16*2, height = 9*2, unit = "cm")


pp2<-(p5|p6)/(p7|p8) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
pp2


ggsave(pp2, filename = "Figures/Properties/FigureB.2.pdf", width = 16*2, height = 9*2, unit = "cm")
