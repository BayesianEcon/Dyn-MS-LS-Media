
rm(list = ls())

#load the libraries

list.of.packages <- c("mvtnorm", "Rcpp", "RcppDist", "RcppParallel", "RcppArmadillo","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(mvtnorm)
library(Rcpp)
library(RcppDist)
library(RcppParallel)
library(RcppArmadillo)
library(data.table)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/RepositoryJASA/")
#####################################

########CHANGE YOUR PATH ###########
setwd("~/Desktop/RepositoryJASA/")
#####################################
### expected running time approximately 45 mins on a MacBook Air M2


sourceCpp("Code/05-Static/Rcpp_rf_it.cpp")
load("Data/Static/Data_Env_single_IT.RData")

unique(DBplane$fb_name)

pw_names<-c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" )

# 
DBplane$National_1<- c(0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,1,1,1,1)
DBplane$National_2<- c(0,1,1,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,1,1,1,1)



EL_princ<-EL_princ[order(EL_princ$ith,EL_princ$jth),]

# library(plyr)

#Main Dimensions
N<- length(EL_x$ith) #number of edges
M<- length(EL_x$jth)/2 #number of unique edges
Npages<- length(DBplane$i) #number of unique pages
K<-2

DBplane<-DBplane[order(DBplane$i),]

# mod1<-lm(DBplane$leaning~DBplane$National_1)
# DBplane$leaning<-mod1$residuals

DBplane$leaning<- (DBplane$leaning - min(DBplane$leaning -0.01))/(max(DBplane$leaning + 0.01) - min(DBplane$leaning -0.0001))
# DBplane$t<- as.numeric(as.factor(DBplane$date))

############## SETUP #######################

pages_names<- c(unique(EL_x$i)) # a vector of page names
beta  = rep(0, Npages)


#rm(EL)

#Hyperparameters for Normal-IW theta

#parameter for Beta 
#gamma prior
a_phi <- 0.001  #1400*100  #1400*150  
b_phi <- 0.001  #0.095

#normal prior
a_gamma_0 <-  -0.1
b_gamma_0 <-   15

#normal prior
a_gamma_1 <- 0.1
b_gamma_1 <- 15

#parameter for poisson 
#gamma prior
a_tau<-  0.001
b_tau <- 0.001

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

zi_a[,1] <- rnorm(Npages, 0, 0.001)
zi_a[,2] <- rnorm(Npages, 0, 0.001)

zi_a_it<-list()


#initialization phi

phi <-  120 #200
#initialization gamma_0

gamma_0 <-  -1 #-0.5

#initialization gamma_1

gamma_1 <-  0.5 #0.3

#initialization tau

tau<- 110

phi_gamma_tau_ite<-list()

zi_a_ref<- cbind(rnorm(Npages, 0,1), rnorm(Npages, 0,1))
zi_b_ref<- cbind(rnorm(Npages, 0,1), rnorm(Npages, 0,1))

zi_a_ref_c<- apply(zi_a_ref,2, FUN = function(x){x - mean(x)})
zi_a_ref_c_mean <- apply(zi_a_ref,2, FUN = function(x){mean(x)})

zi_b_ref_c<- apply(zi_b_ref,2, FUN = function(x){x - mean(x)})
zi_b_ref_c_mean <- apply(zi_b_ref,2, FUN = function(x){mean(x)})


rownames(EL_princ)<-1:nrow(EL_princ)

# 
# lam_ad = rep(2.8 , M)
# mu_mat = matrix(0, M, 3)
# 
# id<-diag(c(0.001,0.001, 0.001))
# Sigma_ad <- rep(list(id),M)


lam_ad_za = rep(  1, Npages)
mu_mat_za = matrix(0, Npages, 2)

id<-diag(c(0.1,0.1))
Sigma_ad_za <- rep(list(id),Npages)

##

lam_ad_zb = rep( 2.8, Npages)
mu_mat_zb = matrix(0, Npages, 2)

id<-diag(c(0.001,0.001))
Sigma_ad_zb <- rep(list(id),Npages)
# 
# DBplane<-DBplane[order(DBplane$t, DBplane$i ), ]
# EL_princ<-EL_princ[order(EL_princ$t, EL_princ$ith ), ]
# EL_x<-EL_x[order(EL_x$t, EL_x$ith ), ]


############

start_time <- Sys.time()

result<- MCMC(
  princ_w = EL_princ$w,
  princ_ones = EL_princ$ones,
  x_w = EL_x$w,
  x_ones =  EL_x$ones,
  mu_ad = c(0,0),
  lm = 1,
  sigma_ad = 1*diag(2),
  leaning = DBplane$leaning, 
  ncomments = DBplane$ncomments ,
  lam_ad_za = lam_ad_za,
  mu_mat_za = mu_mat_za,
  Sigma_ad_za = Sigma_ad_za, 
  lam_ad_zb = lam_ad_za,
  mu_mat_zb = mu_mat_za,
  Sigma_ad_zb = Sigma_ad_za, 
  alpha = 0,
  beta  = rep(0, Npages),
  xi_state1 =  c(1), 
  xi_state2 =  c(0),
  zi_a_1 = zi_a$X1, 
  zi_a_2 = zi_a$X2, 
  zi_b_1 = zi_a$X1, 
  zi_b_2 = zi_a$X2, 
  mu_alpha =0,
  mu_beta =0,
  sigma_alpha = 15,
  sigma_beta = 15,
  mu_a = mu_a , 
  sigma_a = sigma_a,  
  mu_b = mu_a ,  
  sigma_b = sigma_a,
  phi = phi,
  gamma_0 = gamma_0, 
  gamma_1 = gamma_1,
  a_phi = a_phi, 
  b_phi = b_phi, 
  a_gamma_0 = a_gamma_0,
  b_gamma_0 = b_gamma_0, 
  a_gamma_1 = a_gamma_1, 
  b_gamma_1 = b_gamma_1, 
  omega_lower_a = 2, 
  omega_lower_b = 2, 
  P = matrix(0,2,2) ,
  M = M,
  N = Npages,
  Time = 1 ,
  x_i = EL_x$ith,
  DBplane_i  = DBplane$i,
  Iterations = 15000)
end_time <- Sys.time()
end_time-start_time

save(result, file = "Code/05-Static/RESULT_single_it.RData")

##################
######PLOTS#######


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
setwd("~/Desktop/RepositoryJASA/")
#####################################


load("Code/05-Static/RESULT_single_it.RData")
load("Data/Static/Data_Env_single_IT.RData")
source("Code/05-Static/Misc/CoordCartesian.R")

pages_names_it<- c(unique(EL_x$i)) # a vector of page names

top_names<-c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" )

beta_it<- result[[1]]

zeta_it<- result[[3]]
zeta_it<-data.frame(zeta_it)
colnames(zeta_it) <-c("X1","X2", "i", "it")

zeta_it_sub <-zeta_it[zeta_it$it %in% seq(2000, 15000, 1), ]
zeta_it<-zeta_it[,-2]
# check<-reshape(zeta_it, idvar = "it", timevar = "i", direction = "wide")
# check<-check[,-1]
# plot(c(as.matrix(check[2000:15000,])), c(beta_it[2000:15000,]), pch = ".")

agg_zeta_it<-aggregate(zeta_it_sub[,1:2], by = list(zeta_it_sub$i), FUN = mean )
agg_beta_it<-colMeans(beta_it[seq(2000, 15000, 1), ])
agg_zeta_it$names<- pages_names_it



db_it<-data.frame(agg_zeta_it= agg_zeta_it$X1,agg_beta_it= agg_beta_it, names=pages_names_it )
p3<-ggplot(db_it, aes(x = agg_zeta_it, y = agg_beta_it, label = names))
p3<- p3+ geom_point(col = "#ff420f", size = 1) 
p3<- p3 + geom_text_repel(data = db_it[db_it$names %in% top_names, ] ,size = 2,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect" ,  title = "Italy", col = "white")
p3 <- p3 + theme(panel.grid = element_blank())  + xlim(-2.5,2.5) +ylim(-2,7)
p3 <- p3 + theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
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



p3

ggsave(p3, filename = "Figures/Static/latent_static_fe_it.pdf", units = "cm", width = 16*2, height = 9*2 )


result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)
result3$gamma_1<-(result3$gamma_1)


gat1<-gather(result3[result3$ite %in% seq(5000, 15000, 5 ),1:4], key = "parameter", value = "value")


gat<-gat1

gat$parameter<-factor(gat$parameter)
gat$parameter<-factor(gat$parameter, labels = c("gamma[0]", "gamma[1]", "phi", "tau"))

gat<-gat[!gat$parameter %in% c("beta", "tau"), ]

make_label <- function(value) {
  x <- as.character(value)
  bquote(italic(.(x))~subjects)
}

plot_labeller <- function(variable, value) {
  do.call(expression, lapply(levels(value), make_label))
}



colnames(gat)[1]<-"level"

z <- ggplot(gat,  aes(x=value) )
z<- z + geom_histogram( aes(y = after_stat(density)),fill = "#ff420f", alpha = 0.6 , col ="black", bins = 20)
z<- z + geom_density(linetype = 1, adjust = 1.5, n = 10000)
z<-z+ facet_wrap( .~ level, ncol = 6, scales = "free",labeller = label_parsed)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "phi") , fun = dgamma, args = list(0.01, 0.01), geom = "area", fill ="black", alpha = 0.3, linetype = "dashed", colour = "black")
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[0]"),fun = dnorm, args = list(mean = 0, sd = 5), geom = "area",fill ="black", alpha = 0.3, linetype = "dashed", colour = "black")
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[1]"),fun = dnorm, args = list(mean = 0, sd = 5), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black")
z<-z+geom_hline(yintercept = -0.001, color ="white", size = 1.1)
z<- z + scale_shape_manual(labels =  parse_format())
z <- z + labs(title = "Italy", 
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




w_it<-z +  coord_cartesian_panels(
  panel_limits = tibble::tribble(
    ~level, ~xmin, ~xmax,
    "gamma[0]"     ,   mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]),
    "gamma[1]"     , mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]),
    "phi"     ,   0,     mean(gat$value[gat$parameter == "phi" & gat$type == "posterior" ]) + 7*sd(gat$value[gat$parameter == "phi" & gat$type == "posterior" ])
    
  ))

w_it

pew_it<- agg_zeta_it[agg_zeta_it$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ),]
pew_it$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5)
pew_it$score_b<-c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)

#pew_it$X1<- pew_it$X1 -  mean(pew_it$X1)
pew_it$score_a<- pew_it$score_a -  mean(pew_it$score_a)



