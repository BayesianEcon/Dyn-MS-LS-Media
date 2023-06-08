############## Gibbs Sampler #######################
############## Preliminaries #######################


rm(list = ls())

#load the libraries

list.of.packages <-
  c("mvtnorm", "Rcpp", "RcppDist", "RcppParallel", "RcppArmadillo")
new.packages <-
  list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)


library(mvtnorm)
library(Rcpp)
library(RcppDist)
library(RcppParallel)
library(RcppArmadillo)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

Multi <- function(x) {
  m <- rmultinom(1, size = 1, prob = x)
  m <- t(m)
  return(m)
}

odd <- function(x) {
  (1 / (1 + exp(-1 * x)))
} #logistic function
odd_inv <- function(x) {
  (1 - (1 / (1 + exp(-1 * x))))
} #

log_sum_exp = function(x) {
  b <-
    max(x)
  return(b + log(sum(exp(x - b))))
} #define log sum exp function
softmax <-
  function (x) {
    exp(x - log_sum_exp(x))
  } #define softmax function


load("Data/Dynamic/DataEnv_DE_all.RData")
sourceCpp("Code/Model/MS_LS_FE.cpp")


a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )


for(i in 1:length(unique(EL_x$i) )){
  
  name = unique(EL_x$i)[i]
  resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
  x<- rle(resa$x==0)
  
  if(length(x$lengths[x$values == TRUE])>0){
    a[i,2]= max(x$lengths[x$values == TRUE])}
  print(i)
}


thr <- c(a$name[a$max_czeros > 15])

DBplane<-DBplane[! DBplane$fb_name %in% thr ,  ]
DBplane$i<- as.numeric(factor(DBplane$fb_name))
#
EL_x<-EL_x[! EL_x$i %in% thr ,  ]
EL_x<-EL_x[! EL_x$j %in% thr ,  ]
#
EL_princ<-EL_princ[! EL_princ$i %in% thr ,  ]
EL_princ<-EL_princ[! EL_princ$j %in% thr ,  ]
#
EL<-EL[! EL$i %in% thr ,  ]
EL<-EL[! EL$j %in% thr ,  ]

EL$ith<-as.numeric(factor(EL$i) , levels = unique(EL$i))
EL$jth<-as.numeric(factor(EL$j) , levels = unique(EL$i))

EL_x$ith<-as.numeric(factor(EL_x$i) , levels = unique(EL_x$i))
EL_x$jth<-as.numeric(factor(EL_x$j) , levels = unique(EL_x$i))

EL_princ$ith<-as.numeric(factor(EL_princ$i , levels = unique(EL_x$i)))
EL_princ$jth<-as.numeric(factor(EL_princ$j, levels = unique(EL_x$i)))


###########################
#load("DataEnv_DE.RData")

#source("/home/peruzzi/LPmodel/political_programs/02_DE/InputNodesModification_de.R")


EL_x$ith_jth<-paste0(EL_x$ith,"_",EL_x$jth)
EL_x$i_j<- EL_x$ith_jth
EL_x$t<-as.numeric(as.factor(EL_x$t))
DBplane$t<-as.numeric(as.factor(DBplane$t))
EL_princ$t<-as.numeric(as.factor(EL_princ$t))
DBplane$ith<-as.numeric(as.factor(DBplane$i))
DBplane$leaning<- (DBplane$leaning - (min(DBplane$leaning)-0.01))/((max(DBplane$leaning)+0.01) - (min(DBplane$leaning)-0.01)   )


############## SETUP #######################


pages_names<- c(unique(EL$i)) # a vector of page names

#Main Dimensions
N<- length(unique(EL_x$i_j)) #number of edges
M<- length(unique(EL_x$i_j))/2 #number of unique edges
Time<- length(unique(EL_x$t)) #number of unique dates
Npages<- length(pages_names) #number of unique pages
K <- 2 #number of states

#rm(EL)

#parameter for Beta 
#gamma prior
a_phi <- 0.001  #1400*100  #1400*150  
b_phi <-  0.001  #0.095

#normal prior
a_gamma_0 <-  -0.1
b_gamma_0 <-   10

#normal prior
a_gamma_1 <- 0.1
b_gamma_1 <- 10

#Hyperparameters for zi_ik
mu_a<-c( 0 , 0)
mu_b<-c( 0, 0)

sigma_a<- 15
sigma_b<- 15


#Hyperparameters for P
omega_lower_a <- 2
omega_lower_b <- 2

# #Initialization of Z_ia
zi_a<-matrix(c( 0 , 0), nrow= Npages, ncol = 2, byrow=TRUE)
rownames(zi_a)<- pages_names
zi_a<-data.frame(zi_a)
zi_a$i<-pages_names

zi_a[,1] <- rnorm(length(zi_a[,1]), 0, 0.8)
zi_a[,2] <- rnorm(length(zi_a[,2]), 0, 0.8)

#Initialization of Z_ib
zi_b<-matrix(c(0,0), nrow=Npages, ncol = 2, byrow=TRUE)
zi_b[,2]<-rnorm(Npages,mean= 1, sd = sigma_b )
rownames(zi_b)<- pages_names

zi_b<-data.frame(zi_b)
zi_b$i<-pages_names
# 
zi_b[,1] <-  rnorm(length(zi_b[,1]), 0, 0.2)
zi_b[,2] <-  rnorm(length(zi_b[,2]), 0, 0.2)
zi_b_it<-list()



#Initialization of Xi
# 

xi <- data.frame(matrix(rep(1/2,Time*K), nrow =  Time, ncol = K , byrow=TRUE))
xi<-t(apply(xi,1, FUN = Multi ))

xi<-data.frame(xi)
xi$t<-1:nrow(xi)
#xi<-cbind(xi,as.numeric(unique(EL_princ$t)))
colnames(xi)<-c("state1", "state2", "t")


#Initialization of Pi 
P = matrix(c(0.9, 0.1, 0.1, 0.9), nrow =  2, ncol = K , byrow=TRUE)
Plist_it<-list()


#initialization phi

phi <-  10 #200
#initialization gamma_0

gamma_0 <-  -1 #-0.5

#initialization gamma_1

gamma_1 <-  0.5 #0.3



rownames(EL_princ)<-1:nrow(EL_princ)

######


lam_ad_beta = rep(  0.1, Npages)
mu_mat_beta = matrix(0, Npages, 2)

id<-diag(c(1,1))
Sigma_ad_beta <- rep(list(id),Npages)

# 
# id<-diag(c(1))
# Sigma_ad_beta <- rep(list(id),Npages)

########

lam_ad_za = rep(  0.1, Npages)
mu_mat_za = matrix(0, Npages, 2)

id<-diag(c(1,1))
Sigma_ad_za <- rep(list(id),Npages)

lam_ad_zb = rep( 0.1, Npages)
mu_mat_zb = matrix(0, Npages, 2)

id<-diag(c(1, 1))
Sigma_ad_zb <- rep(list(id),Npages)



rg_eq = 0
interp_eq = 1
ms_eq = 1


#################################
start_time <- Sys.time()  

result = MCMC(princ_w = EL_princ$w,
              princ_ones = EL_princ$ones,
              x_w = EL_x$w,
              x_ones =  EL_x$ones,
              leaning = DBplane$leaning, 
              lam_ad_beta = lam_ad_beta,
              mu_mat_beta = mu_mat_beta,
              Sigma_ad_beta = Sigma_ad_beta, 
              lam_ad_za = lam_ad_za,
              mu_mat_za = mu_mat_za,
              Sigma_ad_za = Sigma_ad_za, 
              lam_ad_zb = lam_ad_zb,
              mu_mat_zb = mu_mat_zb,
              Sigma_ad_zb = Sigma_ad_zb, 
              beta = rnorm(Npages, 0 , 0.01),
              xi_state1 = xi$state1,
              zi_a_1 = rnorm(Npages, 0 , 0.01), 
              zi_b_1 = rnorm(Npages, 0 , 0.01), 
              mu_beta =0,
              sigma_beta = 15,
              mu_a = 0 , 
              sigma_a = 15,  
              mu_b = 0 ,  
              sigma_b = 15,
              s_a = 0.1,
              s_b = 0.1,
              phi = 20,
              gamma_0 = 0, 
              gamma_1 = 0,
              a_phi = a_phi, 
              b_phi = b_phi, 
              a_gamma_0 = 0,
              b_gamma_0 = 15, 
              a_gamma_1 = 0, 
              b_gamma_1 = 15, 
              omega_lower_a = 2, 
              omega_lower_b = 2, 
              P = P,
              N = Npages,
              Time = Time ,
              x_i = EL_x$ith,
              DBplane_i  = DBplane$i,
              prop_sd_gamma = 0.007,
              prop_sd_phi = 2, 
              acc_beta = 0.25, 
              acc_zeta_a = 0.25 , 
              acc_zeta_b = 0.25 , 
              pivot =  which(pages_names == "Bild"),
              sign = 1,
              interp_eq = interp_eq, 
              ms_eq = ms_eq,
              rg_eq = rg_eq,
              Iterations = 35000)
end_time <- Sys.time()
end_time-start_time

if(interp_eq ==1 & ms_eq == 1){
  
  save(result, file = "Code/06-Dynamic/RESULT_SC_DEz_fe.RData")
  
  
}else if(interp_eq ==1 & ms_eq == 0){
  
  save(result, file = "Code/06-Dynamic/RESULT_SC_DEz_fe_no_ms.RData")
  
}else if(interp_eq ==0 & ms_eq == 1){
  
  save(result, file = "Code/06-Dynamic/RESULT_SC_DEz_fe_no_interp.RData")
  
}


################Plots########################

list.of.packages <- c( "ggplot2", "ggrepel", "tidyr", "dplyr", "data.table", "patchwork", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")

library(ggplot2)
library(ggrepel)
library(tidyr)
library(scales)
library(dplyr)
library(data.table)
library(patchwork)
library(ggpubr)


########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################


load("Data/Dynamic/DataEnv_DE_all.RData")
load("Code/06-Dynamic/Results/RESULT_SC_DEz_fe.RData")

a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )

for(i in 1:length(unique(EL_x$i) )){
  
  name = unique(EL_x$i)[i]
  resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
  x<- rle(resa$x==0)
  
  if(length(x$lengths[x$values == TRUE])>0){
    a[i,2]= max(x$lengths[x$values == TRUE])}
  print(i)
}


thr <- c(a$name[a$max_czeros > 15])




DBplane<-DBplane[! DBplane$fb_name %in% thr ,  ]
DBplane$i<- as.numeric(factor(DBplane$fb_name))
#
EL_x<-EL_x[! EL_x$i %in% thr ,  ]
EL_x<-EL_x[! EL_x$j %in% thr ,  ]
#
EL_princ<-EL_princ[! EL_princ$i %in% thr ,  ]
EL_princ<-EL_princ[! EL_princ$j %in% thr ,  ]
#
EL<-EL[! EL$i %in% thr ,  ]
EL<-EL[! EL$j %in% thr ,  ]

EL$ith<-as.numeric(factor(EL$i) , levels = unique(EL$i))
EL$jth<-as.numeric(factor(EL$j) , levels = unique(EL$i))

EL_x$ith<-as.numeric(factor(EL_x$i) , levels = unique(EL_x$i))
EL_x$jth<-as.numeric(factor(EL_x$j) , levels = unique(EL_x$i))

EL_princ$ith<-as.numeric(factor(EL_princ$i , levels = unique(EL_x$i)))
EL_princ$jth<-as.numeric(factor(EL_princ$j, levels = unique(EL_x$i)))


dates<- unique(as_date(DBplane$t))


pages_names_de<- c(unique(EL_x$i)) # a vector of page names
pages_names_de[10]<- "FAZ"

beta_de<- result[[1]]

zeta_de_a<- result[[3]]
zeta_de_a<-data.frame(zeta_de_a)

zeta_de_b<- result[[4]]
zeta_de_b<-data.frame(zeta_de_b)


colnames(zeta_de_a) <-c("X1","X2", "i", "it")
colnames(zeta_de_b) <-c("X1","X2", "i", "it")


zeta_de_a_sub <-zeta_de_a[zeta_de_a$it %in% seq(20000, 35000, 1), ]
zeta_de_b_sub <-zeta_de_b[zeta_de_b$it %in% seq(20000, 35000, 1), ]



zi_a_it<-data.frame(result$zi_a_it)
colnames(zi_a_it)<-c("value", "value2", "i", "it")
zi_a_it<-zi_a_it[,-2]

zi_a_it<-reshape2::dcast(zi_a_it, it ~ i)
zi_a_it<-zi_a_it[,-1]

zi_b_it<-data.frame(result$zi_b_it)
colnames(zi_b_it)<-c("value", "value2", "i", "it")
zi_b_it<-zi_b_it[,-2]

zi_b_it<-reshape2::dcast(zi_b_it, it ~ i)
zi_b_it<-zi_b_it[,-1]


# obj1<-update_zi(as.matrix(zi_a_it), as.matrix(zi_b_it))
# zi_a_it<-obj1[[1]]
# zi_b_it<-obj1[[2]]


agg_zeta_de_a<-aggregate(zeta_de_a_sub[,1:2], by = list(zeta_de_a_sub$i), FUN = mean )
agg_zeta_de_a$X1<-colMeans(zi_a_it[seq(20000, 35000, ), ])

agg_beta_de<-colMeans(beta_de[seq(20000, 35000, 1), ])
agg_zeta_de_a$names<- pages_names_de

agg_zeta_de_b<-aggregate(zeta_de_b_sub[,1:2], by = list(zeta_de_b_sub$i), FUN = mean )
agg_zeta_de_b$X1<-colMeans(zi_b_it[seq(20000, 35000, ), ])

agg_beta_de<-colMeans(beta_de[seq(20000, 35000, 1), ])
agg_zeta_de_b$names<- pages_names_de

top_names<-  c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" )

pew_de_a<- agg_zeta_de_a[agg_zeta_de_a$names %in% c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
pew_de_a$State<- "State H"
pew_de_a$score_a<- c(3.6, 3.2, 3.2, 2.7, 3) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
pew_de_a$score_b<- c(3.1, 2.9, 3, 2.8, 2.8) - mean(c(3.1, 2.9, 3, 2.8, 2.8))

cor(pew_de_a$X1, pew_de_a$score_a)

pew_de_b<- agg_zeta_de_b[agg_zeta_de_b$names %in% c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
pew_de_b$State<- "State L"
pew_de_b$score_a<- c(3.6, 3.2, 3.2, 2.7, 3) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
pew_de_b$score_b<- c(3.1, 2.9, 3, 2.8, 2.8) - mean(c(3.1, 2.9, 3, 2.8, 2.8))

pew_de<-rbind(pew_de_a, pew_de_b)
cor(pew_de_b$X1, pew_de_b$score_a)
pew_de$country<- "Germany"

db_de_a<-data.frame(agg_zeta_de= agg_zeta_de_a$X1,agg_beta_de= agg_beta_de, names=pages_names_de )
db_de_a$State<- "State H"
db_de_b<-data.frame(agg_zeta_de= agg_zeta_de_b$X1,agg_beta_de= agg_beta_de, names=pages_names_de )
db_de_b$State<- "State L"



db_de<- rbind(db_de_a,db_de_b)
db_de$State<-factor(db_de$State,  levels = c("State L", "State H"))


p2<-ggplot(db_de, aes(x = agg_zeta_de, y = agg_beta_de, label = names)) +facet_wrap(~State, scales = "free_y")
p2<- p2+ geom_point(col = "#ff420f", size = 1)
p2<- p2 + geom_text_repel(data = db_de,size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Germany", col = "white")
#p2<- p2 + geom_text_repel(data = db_de[db_de$names %in% top_names,] ,size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Germany", col = "white")
p2 <- p2 + theme(panel.grid = element_blank())
p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
p2 <- p2 + theme_minimal() + theme_minimal() + theme(strip.placement = "outside",legend.position="none", 
                                                     plot.title = element_text(face = "italic",   size = 14),
                                                     strip.text.x = element_blank(),
                                                     
                                                     axis.title = element_text( face = "italic",size = rel(1)),
                                                     axis.title.y = element_text(size = 14, angle=90,vjust =2),
                                                     axis.title.x = element_blank(),
                                                     axis.text=element_text(size=12),
                                                     axis.line = element_blank(),
                                                     axis.ticks = element_line(),
                                                     panel.grid = element_blank(),
                                                     panel.border = element_rect(colour = "black", fill = NA))



p2_de<-p2





result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)


gat1<-gather(result3[result3$ite %in% seq(20000, 35000, 10 ),1:4], key = "parameter", value = "value")

result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)


gat2<-gather(result8[result8$ite %in% seq(30000, 50000, 1 ),1:2], key = "parameter", value = "value")
gat2$parameter<-factor(gat2$parameter, levels=c("X2","X1"), labels = c("sigma[L]^2","sigma[H]^2"))

gat1$parameter<-factor(gat1$parameter)
gat1$parameter<-factor(gat1$parameter, labels = c("gamma[0]", "gamma[1]", "phi", "tau"))

gat<-rbind(gat1, gat2)
gat<-gat[!gat$parameter %in% c("alpha", "beta" , "tau"), ]



make_label <- function(value) {
  x <- as.character(value)
  bquote(italic(.(x))~subjects)
}

plot_labeller <- function(variable, value) {
  do.call(expression, lapply(levels(value), make_label))
}


gat<-gat[,1:2]

# gat_prior<-gat
# gat_prior$value <-c(prior_phi, prior_g0, prior_g1)
# gat_prior$type <- "prior"


colnames(gat)[1]<-"level"

gat$type <- "posterior"
gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[H]^2"))


sturges <- function(x){ pretty(range(x),
                               n = nclass.Sturges(x),
                               min.n = 1)}

z <- ggplot(gat,  aes(x=value ))
z <-  z+ geom_histogram(data = gat[gat$level == "gamma[0]",], aes( y = after_stat(density)), binwidth=0.1, fill = "#ff420f", alpha = 0.6 , col ="black", breaks = sturges(gat[gat$level == "gamma[0]","value"] ) ) 
z <- z+ geom_histogram(data = gat[gat$level == "gamma[1]",],  aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "gamma[1]","value"] ) ) 
z <- z+ geom_histogram(data = gat[gat$level == "phi",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "phi","value"] ) ) 
z <- z+ geom_histogram(data = gat[gat$level == "sigma[L]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[L]^2","value"] ) ) 
z <- z+ geom_histogram(data = gat[gat$level == "sigma[H]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[H]^2","value"] ) ) 
z<- z + geom_density(linetype = 1,  n = 10000, adjust = 2)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "phi") , fun = dgamma, args = list(shape=0.01, rate =0.01), geom = "area", fill ="black", alpha = 0.3, linetype = "dashed", colour = "black", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[0]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3, linetype = "dashed", colour = "black", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[1]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[L]^2"),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[H]^2"),fun = dinvgamma, args = list(shape = 0.01, scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed) 
z<- z + scale_shape_manual(labels =  parse_format())
z <- z + labs(title = "Germany", 
              x = "Value",  y = "")+  theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                                                size = 22, hjust = 0.5),
                                                              axis.text=element_text(size=8),
                                                              strip.text.x = element_text(size = 18, face = "italic"),
                                                              axis.title.x = element_blank(),
                                                              axis.title.y = element_text(size = 22, angle=90,vjust =2, face = "italic"),
                                                              axis.line = element_blank(),
                                                              axis.ticks = element_line(),
                                                              axis.ticks.x  = element_line(),
                                                              panel.spacing.x = unit(4, "mm"),
                                                              panel.grid = element_blank(),
                                                              panel.border = element_rect(colour = "black", fill = NA))




w_de<-z 


result7<- result[[6]]
colnames(result7)<-c("state1", "state2", "t", "ite")
result7<-data.frame(result7)


xi_it_sub<-result7[result7$ite %in% seq(25000, 35000, 10),]
xi_it_agg<-aggregate(xi_it_sub[,1:2], by = list(t = xi_it_sub$t), FUN = mean)
xi_it_agg$State<-  xi_it_agg$state1
# xi_it_agg$State<- factor(xi_it_agg$state1)
# xi_it_agg$State<- factor(xi_it_agg$State, labels = c("L", "H"))
xi_it_agg$t <- dates

xi_it_agg_area<- xi_it_agg
xi_it_agg_area$dd<- c(1, diff(xi_it_agg_area$state1))
xi_it_agg_area<-xi_it_agg_area[(abs(xi_it_agg_area$dd) >=  0.20 ) ,]
xi_it_agg_area$lead<-lead(xi_it_agg_area$t)
xi_it_agg_area$lead[length(xi_it_agg_area$lead)]<- xi_it_agg$t[length(xi_it_agg$t)]+1

q<- ggplot(xi_it_agg)+labs(x ="time", y = "State", colour = "State")
#q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
q<- q + geom_rect(data = xi_it_agg_area , aes(xmin = t, xmax = lead,  ymax  = state1),  ymin= min(as.numeric(xi_it_agg$state1)), fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
#q<- q + geom_step(aes(x = t, y = as.numeric(State)), color ="black", shape = 15 , size = 0.2)

# q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
#                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# 
q<- q + labs(title = "Germany",
             x = "Time", y = "P.P. - State H")+ theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                                                            size = 36, hjust = 0.5),
                                                                        strip.text.x = element_text(size = 24,face = "italic"),
                                                                        
                                                                        axis.title = element_text( face = "italic",size = rel(1)),
                                                                        axis.title.y = element_text(size = 22, angle=90,vjust =2),
                                                                        axis.title.x = element_text(size = 22, vjust = -0.2),
                                                                        axis.text=element_text(size=16),
                                                                        axis.line = element_blank(),
                                                                        axis.ticks = element_line(),
                                                                        panel.grid = element_blank(),
                                                                        legend.title = element_text(face="italic", size=16),
                                                                        legend.text = element_text(size=14),
                                                                        panel.border = element_rect(colour = "black", fill = NA))



q_de<- q





result6<- data.frame(result[[5]])
result6<- result6[seq(25000, 35000, 10),1:4]
colnames(result6)<-c("p11","p12", "p21","p22")
ps<- apply(result6,2, FUN = mean)
ps
ps_m<-matrix(ps, 2,2, byrow = T)
ps_m<-round(ps_m,3)



melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
melted_cormat$value<-round(melted_cormat$value,3)
melted_cormat$Var1<-factor(melted_cormat$Var1)
melted_cormat$Var2<-factor(melted_cormat$Var2)

melted_cormat$Var1 <- relevel(melted_cormat$Var1, "2")

melted_cormat$Var1 <- factor(melted_cormat$Var1 ,labels = c("L","H"))
melted_cormat$Var2 <- factor(melted_cormat$Var2 ,labels = c("H","L"))


# Create a ggheatmap
ggheatmap1 <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+ labs(x= "States", y ="States", title = "Germany" )+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ coord_flip()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                                                                 size = 22, hjust = 0.5),
                                                                               strip.text.x = element_text(size = 24,face = "italic"),
                                                                               
                                                                               axis.title = element_text( face = "italic",size = rel(1)),
                                                                               axis.title.y = element_text(size = 22, angle=90,vjust =2),
                                                                               axis.title.x = element_text(size = 22, vjust = -0.2),
                                                                               axis.text=element_text(size=16),
                                                                               axis.line = element_blank(),
                                                                               axis.ticks = element_line(),
                                                                               panel.grid = element_blank(),
                                                                               panel.border = element_rect(colour = "black", fill = NA))



ggheatmap_de <- ggheatmap1



