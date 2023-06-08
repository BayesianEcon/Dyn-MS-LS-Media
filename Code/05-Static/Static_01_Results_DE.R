
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
setwd("~/Desktop/Repository/")
#####################################
### expected running time approximately 45.19719 mins on a MacBook Air M2


sourceCpp("Code/05-Static/Rcpp_rf_de.cpp")
load("Data/Static/Data_Env_single_DE.RData")


EL_x$ith_jth<-paste0(EL_x$ith,"_",EL_x$jth)
EL_x$i_j<- EL_x$ith_jth

DBplane$ith<-as.numeric(as.factor(DBplane$i))
DBplane$leaning<- (DBplane$leaning - (min(DBplane$leaning)-0.01))/((max(DBplane$leaning)+0.01) - (min(DBplane$leaning)-0.01)   )


############## SETUP #######################

pages_names<- c(unique(EL_x$i)) # a vector of page names

#Main Dimensions
N<- length(unique(EL_x$i_j)) #number of edges
M<- length(unique(EL_x$i_j))/2 #number of unique edges
Time<- 1 #number of unique dates
Npages<- length(pages_names) #number of unique pages
K <- 2 #number of states

DBplane<-DBplane[order(DBplane$i),]

# mod1<-lm(DBplane$leaning~DBplane$National_1)
# DBplane$leaning<-mod1$residuals

# DBplane$t<- as.numeric(as.factor(DBplane$date))

############## SETUP #######################

pages_names<- c(unique(EL_x$i)) # a vector of page names
beta  = rep(0, Npages)


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

lam_ad_za = rep(  1, Npages)
mu_mat_za = matrix(0, Npages, 2)

id<-diag(c(0.1,0.1))
Sigma_ad_za <- rep(list(id),Npages)

##

lam_ad_zb = rep( 2.8, Npages)
mu_mat_zb = matrix(0, Npages, 2)

id<-diag(c(0.001,0.001))
Sigma_ad_zb <- rep(list(id),Npages)

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

save(result, file = "Code/05-Static/RESULT_single_de.RData")

#####Plots########

#load the libraries

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


load("Code/05-Static/RESULT_single_de.RData")
load("Data/Static/Data_Env_single_DE.RData")
source("Code/05-Static/Misc/CoordCartesian.R")

pages_names_de<- c(unique(EL_x$i)) # a vector of page names

beta_de<- result[[1]]

zeta_de<- result[[3]]
zeta_de<-data.frame(zeta_de)
colnames(zeta_de) <-c("X1","X2", "i", "it")

zeta_de_sub <-zeta_de[zeta_de$it %in% seq(2000, 15000, 1), ]

agg_zeta_de<-aggregate(zeta_de_sub[,1:2], by = list(zeta_de_sub$i), FUN = mean )
agg_beta_de<-colMeans(beta_de[seq(2000, 15000, 1), ])
agg_zeta_de$names<- pages_names_de

zeta_de_sub <-zeta_de[zeta_de$it %in% seq(2000, 15000, 1), ]
zeta_de<-zeta_de[,-2]

top_names<-  c("Bild", "FAZ.NET - Frankfurter Allgemeine Zeitung", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" )

pew_de<- agg_zeta_de[agg_zeta_de$names %in% c("Bild", "FAZ.NET - Frankfurter Allgemeine Zeitung", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
pew_de$score_a<- c(3.6, 3.2, 3.2, 2.7, 3)
pew_de$score_b<- c(3.1, 2.9, 3, 2.8, 2.8)

db_de<-data.frame(agg_zeta_de= agg_zeta_de$X1,agg_beta_de= agg_beta_de, names=pages_names_de )
p1<-ggplot(db_de, aes(x = agg_zeta_de, y = agg_beta_de, label = names))
p1<- p1+ geom_point(col = "#ff420f", size = 1) 
p1<- p1 + geom_text_repel(data = db_de[db_de$names %in% top_names,] ,size = 2,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Germany", col = "white")
p1 <- p1 + theme(panel.grid = element_blank())  + xlim(-2.5,2.5) +ylim(-2,7)
p1 <- p1 + theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                     size = 22, hjust = 0.5),
                                   strip.text.x = element_text(size = 24,face = "italic"),
                                   
                                   axis.title = element_text( face = "italic",size = rel(1)),
                                   axis.title.y = element_text(size = 22, angle=90,vjust =2),
                                   axis.title.x = element_text(size = 22, vjust = -0.2),
                                   axis.text=element_text(size=16),
                                   axis.line = element_blank(),
                                   axis.ticks = element_line(),
                                   axis.ticks.x  = element_line(),
                                   panel.grid = element_blank(),
                                   panel.border = element_rect(colour = "black", fill = NA))



p1


ggsave(p1, filename = "Figures/Static/latent_static_fe_de.pdf", units = "cm", width = 16*2, height = 9*2 )



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


gat<-gat[,1:2]
gat$type <- "posterior"

# gat<-rbind(gat, gat_prior)



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




w_de<-z +  coord_cartesian_panels(
  panel_limits = tibble::tribble(
    ~level, ~xmin, ~xmax,
    "gamma[0]"     ,   mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]),
    "gamma[1]"     , mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]),
    "phi"     ,   0,     mean(gat$value[gat$parameter == "phi" & gat$type == "posterior" ]) + 7*sd(gat$value[gat$parameter == "phi" & gat$type == "posterior" ])
    
  ))

w_de
#pew_de$X1<- pew_de$X1 -  mean(pew_de$X1)
pew_de$score_a<- pew_de$score_a -  mean(pew_de$score_a)



