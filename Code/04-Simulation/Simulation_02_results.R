rm(list = ls())

#load the libraries

list.of.packages <- c("mvtnorm", "Rcpp", "RcppDist", "RcppParallel", "RcppArmadillo","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

########CHANGE YOUR PATH ###########
setwd("~/Documents/GitHub/Dyn-MS-LS-Media/")
#####################################

library(data.table)
library(mvtnorm)
library(Rcpp)
library(RcppDist)
library(RcppParallel)
library(RcppArmadillo)
############## Synthetic Data Generation #######################

source("Code/04-Simulation/Simulation_01_dgp.R")

############## Gibbs Sampler #######################
############## Preliminaries #######################

load("Data/Simulation/SimulationEnv_FE.RData")
sourceCpp("Code/Model/MS_LS_FE.cpp")


Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}


EL_x <- as.data.table(EL[,c("i","w", "wl1","ones", "t", "i_j")])
EL_x$t<-as.numeric(EL_x$t)
setkey(EL_x, t)


DBplane<-data.frame(DBplane)

theta_ij_checklist <- list()
zi_a_checklist     <- list()
zi_b_checklist     <- list()

############## SETUP #######################

xi_1_t <- xi$S1
xi_2_t <- xi$S2


pages_names<- c(unique(EL$i)) # a vector of page names

#Main Dimensions
N<- length(unique(EL_x$i_j)) #number of edges
M<- length(unique(EL_x$i_j))/2 #number of unique edges
Time<- length(unique(EL_x$t)) #number of unique dates
Npages<- length(pages_names) #number of unique pages
K <- 2 #number of states

#rm(EL)

var_a_ref<-var(zi_a_ref[,1])
var_b_ref<-var(zi_b_ref[,1])

#parameter for Beta 
#gamma prior
a_phi <- 0.001  #1400*100  #1400*150  
b_phi <-  0.001  #0.095

#normal prior
a_gamma_0 <-   -0.01
b_gamma_0 <-   5

#normal prior
a_gamma_1 <- 0.01
b_gamma_1 <- 5

#parameter for poisson 
#gamma prior
a_tau<- 0.001
b_tau <- 0.001

#Hyperparameters for zi_ik
mu_a<-c( 0 , 0)
mu_b<-c( 0, 0)

sigma_a<- 10
sigma_b<- 10


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

gamma_0 <-  -0.01 #-0.5

#initialization gamma_1

gamma_1 <-  0.01 #0.3



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

interp_eq = 1 #either 1 interpretation equation or 0
ms_eq = 1
rg_eq = 0

#################################
## expected running time  12.29215 mins on MacBook Air M2

Iterations<- 20000

start_time <- Sys.time()  

result = MCMC(princ_w = EL_princ$w, #edge weights from the network triangle dataset
              princ_ones = EL_princ$ones,  #vector of ones from the network triangle dataset
              x_w = EL_x$w, #edge weights from the network off-diagonal dataset
              x_ones =  EL_x$ones, #vector of ones from the network off-diagonal dataset
              leaning = DBplane$leaning, #text analysis leaning from the dbplane dataset
              lam_ad_beta = lam_ad_beta, #adaptive RW-MH lambda - individual effects
              mu_mat_beta = mu_mat_beta,  #adaptive RW-MH mu - individual effects
              Sigma_ad_beta = Sigma_ad_beta,  #adaptive RW-MH sigma - individual effects
              lam_ad_za = lam_ad_za, #adaptive RW-MH lambda - latent coordinate - state a
              mu_mat_za = mu_mat_za, #adaptive RW-MH mu - latent coordinate - state a
              Sigma_ad_za = Sigma_ad_za,   #adaptive RW-MH sigma - latent coordinate - state a
              lam_ad_zb = lam_ad_zb, #adaptive RW-MH lambda - latent coordinate - state b
              mu_mat_zb = mu_mat_zb, #adaptive RW-MH mu - latent coordinate - state b
              Sigma_ad_zb = Sigma_ad_zb,  #adaptive RW-MH sigma - latent coordinate - state b
              beta = rnorm(Npages, 0 , 0.01), # starting value - individual effect
              xi_state1 = xi[,1],  # starting value - state 1
              zi_a_1 = rep(0, Npages), # starting value - latent coordinate - state a
              zi_b_1 = rep(0, Npages), # starting value - latent coordinate - state b
              mu_beta =0, #prior mean - individual effect
              sigma_beta = 15,  #prior sd - individual effect
              mu_a = 0 ,  #prior mean - latent coordinate - state a
              sigma_a = 10,  #starting sd - latent coordinate - state a
              mu_b = 0 ,   #prior mean - latent coordinate - state b
              sigma_b = 10,  #starting sd - latent coordinate - state b
              s_a = 0.1,  #prior shape parameter of sigma (gamma(s_a, s_b))
              s_b = 0.1,  #prior rate parameter of sigma (gamma(s_a, s_b))
              phi = 50, #starting value phi
              gamma_0 = 0,  #starting value -  gamma_0
              gamma_1 = 0,  #starting value - gamma_1
              a_phi = a_phi, #prior shape parameter of phi (gamma(a_phi, b_phi))
              b_phi = b_phi,  #prior rate parameter of phi (gamma(a_phi, b_phi))
              a_gamma_0 = 0,   #prior mean parameter of gamma_0
              b_gamma_0 = 15,  #prior sd parameter of gamma_0
              a_gamma_1 = 0,    #prior mean parameter of gamma_1
              b_gamma_1 = 15,    #prior sd parameter of gamma_1
              omega_lower_a = 2, #starting value - 1st parameter dirichlet
              omega_lower_b = 2, #starting value - 2nd parameter dirichlet
              P = P, #transition Matrix P
              N = Npages, # Number of nodes
              Time = 100 , #Number of times
              x_i = EL_x$i,  #index vector of nodes in the network off-diagonal dataset
              DBplane_i  = DBplane$i,  #index vector of nodes in the leaning dataset
              prop_sd_gamma = 0.01, #proposal RWMH for the gamma parameters
              prop_sd_phi = 35,  #proposal RWMH for the phi parameter
              acc_beta = 0.25 , # target acceptance Adaptive RW-MH -  individual effects
              acc_zeta_a = 0.25 ,  # target acceptance Adaptive RW-MH -  latent coordinates - state a
              acc_zeta_b = 0.25 , # target acceptance Adaptive RW-MH -  latent coordinates - state b
              pivot = 3,  # pivot node (news outlet) to solve identification
              sign= -1,  # sign of the news outlet -1 left or 1 right
              interp_eq = interp_eq, # 0-1 option for the use of the interpretation equation
              ms_eq = ms_eq, # 0-1 option for the use of the markov switching
              rg_eq = rg_eq, # 0-1 option random graph
              Iterations = Iterations #number of iterations
              )
end_time <- Sys.time()
end_time-start_time


if(rg_eq == 0){

if(interp_eq ==1 & ms_eq== 1){
  
  save(result, file = "Code/04-Simulation/RESULT_Sim.RData")
  
  
}else if(interp_eq ==0 & ms_eq== 1){
  
  save(result, file = "Code/04-Simulation/RESULT_Sim_no_interp.RData")
  
}else if(interp_eq ==1 & ms_eq== 0){
  
  save(result, file = "Code/04-Simulation/RESULT_Sim_no_ms.RData")
  
}else if(interp_eq ==0 & ms_eq== 0){
  
  save(result, file = "Code/04-Simulation/RESULT_Sim_no_ms_no_int.RData")
  
}

  
}else{
  
  save(result, file = "Code/Simulation/RESULT_Sim_rg.RData")
  
}


######################## 
########Plotting####### 
# Plots are designed for our Model (options: interp_eq = 1, ms_eq = 1, rg_eq = 0)

load("Code/04-Simulation/RESULT_Sim.RData")

list.of.packages <- c("ggplot2", "ggpubr", "patchwork", "tidyr", "scales")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidyr)
library(scales)

result1<- result[[1]]
result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)

result4<- result[[3]]
result4<-data.frame(result4)
colnames(result4) <-c("X1","X2", "i", "it")

result5<- result[[4]]
result5<-data.frame(result5)
colnames(result5) <-c("X1","X2", "i", "it")

agg_beta = colMeans(result1[,1:20])

zi_a_it_sub<-result4[result4$it %in%  seq(10000, Iterations, 1), ]
agg_zi_a_it<-aggregate(zi_a_it_sub$X1, by = list(zi_a_it_sub$i), FUN = mean)

zi_b_it_sub<-result5[result5$it %in% seq(10000, Iterations, 1), ]
agg_zi_b_it<-aggregate(zi_b_it_sub$X1, by = list(zi_b_it_sub$i), FUN = mean)


result4$state <- "a"
result5$state <- "b"

#credible ellipse state a

DB_a<- data.frame(X1 = 0, X2 = 0, i = 0, state = "a" )

for(i in 1:Npages){
  
  try<- result4[result4$i == i, ]
  try$X2<- result1[,i]
  try<-try[   seq(10000, Iterations, 1), ]
  
  
  CV = cov(try[,1:2])
  mean= c(mean(try$X1), mean(try$X2) )
  ev <- eigen(CV)
  ev$values
  ev$vectors
  
  max_eig_val = max(ev$values)
  min_eig_val = min(ev$values)
  
  max_eig_v =   ev$vectors[, which.max(ev$values)]
  min_eig_v =   ev$vectors[, which.min(ev$values)]
  
  chisquare_val = sqrt(qchisq(0.99, df = 2))
  
  angle = atan2(max_eig_v[2],max_eig_v[1])
  if(angle < 0){angle = angle + 2*pi}
  phi = angle
  
  a=chisquare_val*sqrt(max_eig_val)
  b=chisquare_val*sqrt(min_eig_val)
  
  theta_grid = seq( 0,2*pi, length.out = 100)
  
  ellipse_x_r  = a*cos( theta_grid )
  ellipse_y_r  = b*sin( theta_grid )
  
  R = rbind(c(cos(phi), sin(phi)), c(-1*sin(phi), cos(phi)))
  
  ellipse_r = rbind(ellipse_x_r,ellipse_y_r)
  r_ellipse = t(ellipse_r) %*% R
  r_ellipse[,1] =   r_ellipse[,1] + mean[1]
  r_ellipse[,2] =   r_ellipse[,2] + mean[2]
  
  db_aux = cbind(X1 = r_ellipse[,1], X2 =  r_ellipse[,2], i = i, state = "a"  )
  
  DB_a<-rbind(DB_a,db_aux )
  print(i)
}

DB_a<-DB_a[-1,]

#credible ellipse state b
DB_b<- data.frame(X1 = 0, X2 = 0, i = 0, state = "b" )

for(i in 1:Npages){
  
  try<- result5[result5$i == i, ] 
  try$X2<- result1[,i]
  try<-try[   seq(10000, Iterations, 1), ]
  
  CV = cov(try[,1:2])
  mean= c(mean(try$X1), mean(try$X2) )
  ev <- eigen(CV)
  ev$values
  ev$vectors
  
  max_eig_val = max(ev$values)
  min_eig_val = min(ev$values)
  
  max_eig_v =   ev$vectors[, which.max(ev$values)]
  min_eig_v =   ev$vectors[, which.min(ev$values)]
  
  chisquare_val = sqrt(qchisq(0.99, df = 2))
  
  angle = atan2(max_eig_v[2],max_eig_v[1])
  if(angle < 0){angle = angle + 2*pi}
  phi = angle
  
  a=chisquare_val*sqrt(max_eig_val)
  b=chisquare_val*sqrt(min_eig_val)
  
  theta_grid = seq( 0,2*pi, length.out = 100)
  
  ellipse_x_r  = a*cos( theta_grid )
  ellipse_y_r  = b*sin( theta_grid )
  
  R = rbind(c(cos(phi), sin(phi)), c(-1*sin(phi), cos(phi)))
  
  ellipse_r = rbind(ellipse_x_r,ellipse_y_r)
  r_ellipse = t(ellipse_r) %*% R
  r_ellipse[,1] =   r_ellipse[,1] + mean[1]
  r_ellipse[,2] =   r_ellipse[,2] + mean[2]
  
  db_aux = cbind(X1 = r_ellipse[,1], X2 =  r_ellipse[,2], i = i, state = "b"  )
  
  DB_b<-rbind(DB_b,db_aux )
  print(i)
}

DB_b<-DB_b[-1,]

DBplot<-rbind(DB_a, DB_b)
DBplot<-data.frame(DBplot)
DBplot$X1<-round(as.numeric(DBplot$X1), 5)
DBplot$X2<-round(as.numeric(DBplot$X2), 5)


agg_zi_a_it$X2<- agg_beta
agg_zi_b_it$X2<- agg_beta


zi_a_ref <- data.frame(zi = zi_a_t, beta =  beta_t )
zi_a_ref$State = "State 1"

zi_b_ref <- data.frame(zi = zi_b_t, beta =  beta_t )
zi_b_ref$State = "State 2"

zi_ref<-rbind(zi_a_ref, zi_b_ref)




sub_res_a =  result4[result4$it %in%  seq(10000, Iterations, 1), ]
sub_res_b =  result5[ result5$it %in%   seq(10000, Iterations, 1), ]
DB_points<- data.frame(rbind( sub_res_a, sub_res_b ))

means<-aggregate(DB_points[,1:2], by = list(DB_points$i, DB_points$state), FUN = mean )
means$State<- means$Group.2

means$State<-factor(means$State, labels = c("State H", "State L"))
DBplot$State<-factor(DBplot$state, labels = c("State H", "State L"))
zi_ref$State <-factor(zi_ref$State , labels = c("State H", "State L"))

means$State<-factor(means$State, levels = c("State L", "State H"))
means$beta<-agg_beta
DBplot$State<-factor(DBplot$State, levels = c("State L", "State H"))
zi_ref$State <-factor(zi_ref$State , levels = c("State L", "State H"))

# DBplot$X1_no_int<- DBplot_no_interp$X1
# DBplot$X2_no_int<- DBplot_no_interp$X2


p<-ggplot(DBplot, aes(x= X1 , y= X2 ) ) +
  facet_grid(cols = vars(State))+
  #geom_hline(yintercept = 0, linetype = "dashed") +  geom_vline(xintercept = 0,  linetype = "dashed")+ 
  #geom_polygon( aes(x= X1_no_int , y= X2_no_int , group = i) , color = "black" , linetype = "dashed",  alpha = 0.1, size = 0.1  )  +
  geom_polygon(fill = "#ff420f",aes(group = i), alpha = 0.5, color = "#ff420f", linewidth = 0.01)+
  #geom_point(means, mapping =  aes(x= X1, y =X2, colour = State) , shape = 19,  size = 0.01 )+
  #geom_text_repel(data= means, aes( x= X1 , y= beta, label = Group.1))+
  geom_point(zi_ref, mapping =  aes(x= zi, y =beta), color = "black" , inherit.aes = FALSE , shape = 4,  size = 1  )  +
  
  labs(#title = "Latent Space", 
    x = "Political Leaning", y = "Individual Effect") + theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                                                                    size = 36, hjust = 0.5),
                                                                                strip.text.x = element_text(size = 14,face = "italic"),
                                                                                
                                                                                axis.title = element_text( face = "italic",size = rel(1)),
                                                                                axis.title.y = element_text(size = 14, angle=90,vjust =2),
                                                                                axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                axis.text=element_text(size=16),
                                                                                axis.line = element_blank(),
                                                                                axis.ticks = element_line(),
                                                                                panel.grid = element_blank(),
                                                                                legend.title = element_text(face="italic", size=16),
                                                                                legend.text = element_text(size=14),
                                                                                panel.border = element_rect(colour = "black", fill = NA))


p


result7<- result[[6]]
colnames(result7)<-c("state1", "state2", "t", "ite")
result7<-data.frame(result7)

#load("~/Desktop/political_programs/political_programs/SimulationEnv.RData")

xi_ref <- xi

library(ggplot2)
xi_it_sub<-result7[result7$ite %in% seq(10000, Iterations, 1),]
xi_it_agg<-aggregate(xi_it_sub[,1:2], by = list(t = xi_it_sub$t), FUN = mean)
xi_it_agg$State<- factor(xi_it_agg$state1)
xi_it_agg$State<- factor(xi_it_agg$State, labels = c("L", "H"))


q<- ggplot(xi_it_agg)+labs(x ="time", y = "State", colour = "State")
q<- q + geom_step(aes(x = t, y = State, group = 1), color = "black", linetype = "dashed", inherit.aes = F) +geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 )
q<- q + labs(#title = "Estimated Latent States",
  x = "Time", y = "States")+ theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                                         size = 36, hjust = 0.5),
                                                     strip.text.x = element_text(size = 24,face = "italic"),
                                                     
                                                     axis.title = element_text( face = "italic",size = rel(1)),
                                                     axis.title.y = element_text(size = 14, angle=90,vjust =2),
                                                     axis.title.x = element_text(size = 14, vjust = -0.2),
                                                     axis.text=element_text(size=16),
                                                     axis.line = element_blank(),
                                                     axis.ticks = element_line(),
                                                     panel.grid = element_blank(),
                                                     legend.title = element_text(face="italic", size=16),
                                                     legend.text = element_text(size=14),
                                                     panel.border = element_rect(colour = "black", fill = NA))



q


mc4<- (p|q) + plot_annotation(tag_levels = 'A')

result8<-result$Sigma_z_ite

gat1<-gather(result3[result3$ite %in% seq(10000, Iterations, 1 ),1:4], key = "parameter", value = "value")


result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)


gat2<-gather(result8[result8$ite %in% seq(10000, Iterations, 1 ),1:2], key = "parameter", value = "value")
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

datax<-data.frame(xx= c(0,0,0, 0, 0) , yy =  c(-0.1, 0.5, 200, var(zi_b_ref[,1]),var(zi_a_ref[,1]) )) 
datax$parameter = unique(gat$parameter)[order(unique(gat$parameter))]

library(scales)
z <- ggplot(gat,  aes(y=value)) +
  geom_boxplot(fill = "#ff420f", alpha = 0.7 ) + facet_wrap( ~ parameter, ncol = 6, scales = "free_y",labeller = label_parsed)
z<- z + scale_shape_manual(labels =  parse_format())
z<- z + geom_hline(data= datax ,  aes(yintercept = yy, group = parameter ), pch = 4, linetype = "dashed",  col = "black", inherit.aes = F)
z <- z + labs(#title = "Boxplots Marginal Posteriors", 
  x = "",  y = "Value")+  theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                                    size = 36, hjust = 0.5),
                                                  axis.text=element_text(size=16),
                                                  strip.text.x = element_text(size = 18, face = "italic"),
                                                  axis.title.x = element_blank(),
                                                  axis.text.x = element_blank(),
                                                  axis.title.y = element_text(size = 14, angle=90,vjust =2, face = "italic"),
                                                  axis.line = element_blank(),
                                                  axis.ticks = element_line(),
                                                  axis.ticks.x = element_blank(),
                                                  panel.grid = element_blank(),
                                                  panel.border = element_rect(colour = "black", fill = NA))



z

design<-"11111
         11111
         22222
         22222
         33333"

mc5<- p +q+z + plot_annotation(tag_levels = 'A') + plot_layout(design = design)  &theme(axis.text= element_text(size = 12))
ggsave(mc5, filename = "Figures/Simulation/Figure5.pdf", units = "cm", width = 21, height = 29.7 )



########### Trace Plots ###################################

result1<- data.frame(result[[1]])
result1$X3<- result1$X3+1
result1$mean1 <- cumsum(result1$X1)/(1:length(result1$X1))
result1$mean2 <- cumsum(result1$X2)/(1:length(result1$X2))
result1<-data.frame(result1)
result1$ite<-1:nrow(result1)
result1$X1_mean <- cumsum(result1$X1)/(1:length(result1$X1))
result1$X11_mean <- cumsum(result1$X11)/(1:length(result1$X11))


result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)
result3$phi_mean <- cumsum(result3$phi)/(1:length(result3$phi))
result3$gamma0_mean <- cumsum(result3$gamma_0)/(1:length(result3$gamma_0))
result3$gamma1_mean <- cumsum(abs(result3$gamma_1))/(1:length(result3$gamma_1))
result3$tau_mean <- cumsum(result3$tau)/(1:length(result3$tau))


result4<- result[[3]]
result5<- result[[4]]





p10<-ggplot(result1) + geom_line(aes(x = ite , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = X1_mean), col = "black")  + geom_hline(yintercept  = beta_t[1], linetype = "dashed", color = "black")
p10<- p10 +  labs(title = expression(paste("Parameter ", alpha[1])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                               size = rel(1.2), hjust = 0.5),
                                                                                                                           
                                                                                                                           axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                           axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                           axis.title.x = element_text(vjust = -0.2),
                                                                                                                           axis.line = element_line(colour="black"),
                                                                                                                           axis.ticks = element_line(),
                                                                                                                           legend.title = element_text(face="italic"),
                                                                                                                           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p10


p11<-ggplot(result1) + geom_line(aes(x = ite , y = X11), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = X11_mean), col = "black")  + geom_hline(yintercept  = beta_t[11], linetype = "dashed", color = "black")
p11<- p11 +  labs(title = expression(paste("Parameter ", alpha[11])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                size = rel(1.2), hjust = 0.5),
                                                                                                                            
                                                                                                                            axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                            axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                            axis.title.x = element_text(vjust = -0.2),
                                                                                                                            axis.line = element_line(colour="black"),
                                                                                                                            axis.ticks = element_line(),
                                                                                                                            legend.title = element_text(face="italic"),
                                                                                                                            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p11




p7<-ggplot(result3) + geom_line(aes(x = ite , y = phi), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = phi_mean), col = "black")  + geom_hline(yintercept  = 200, linetype = "dashed", color = "black")
p7<- p7 +  labs(title = expression(paste("Parameter ", phi)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                        size = rel(1.2), hjust = 0.5),
                                                                                                                    
                                                                                                                    axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                    axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                    axis.title.x = element_text(vjust = -0.2),
                                                                                                                    axis.line = element_line(colour="black"),
                                                                                                                    axis.ticks = element_line(),
                                                                                                                    legend.title = element_text(face="italic"),
                                                                                                                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p7


p8<-ggplot(result3) + geom_line(aes(x = ite , y = gamma_0), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = gamma0_mean), col = "black")  + geom_hline(yintercept  = -0.1, linetype = "dashed", color = "black") + ylim( -0.15, 0)
p8<- p8 +  labs(title = expression(paste("Parameter ", gamma[0])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                             size = rel(1.2), hjust = 0.5),
                                                                                                                         
                                                                                                                         axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                         axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                         axis.title.x = element_text(vjust = -0.2),
                                                                                                                         axis.line = element_line(colour="black"),
                                                                                                                         axis.ticks = element_line(),
                                                                                                                         legend.title = element_text(face="italic"),
                                                                                                                         strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p8


p9<-ggplot(result3) + geom_line(aes(x = ite , y =  abs(gamma_1) ), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = gamma1_mean), col = "black")  + geom_hline(yintercept  =  0.5, linetype = "dashed", color = "black")
p9<- p9 +  labs(title = expression(paste("Parameter ", gamma[1])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                             size = rel(1.2), hjust = 0.5),
                                                                                                                         
                                                                                                                         axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                         axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                         axis.title.x = element_text(vjust = -0.2),
                                                                                                                         axis.line = element_line(colour="black"),
                                                                                                                         axis.ticks = element_line(),
                                                                                                                         legend.title = element_text(face="italic"),
                                                                                                                         strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p9


mc2<-(p10|p11)/(p8|p9)/(p7)  + plot_annotation(tag_levels = 'A')
ggsave(mc2, filename = "Figures/Simulation/FigureD1.pdf", units = "cm", width = 16*2, height = 9*2 )



result4<-data.frame(result4)
colnames(result4) <-c("X1","X2", "i", "it")
try<-data.frame(result4)



try<-data.frame(result4)
try = try[try$i == 1, ]
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))


p12<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zi_a_t[1], linetype = "dashed", color = "black")
p12<- p12 +  labs(title = expression(paste("Parameter ", zeta["1, H"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                   size = rel(1.2), hjust = 0.5),
                                                                                                                               
                                                                                                                               axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                               axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                               axis.title.x = element_text(vjust = -0.2),
                                                                                                                               axis.line = element_line(colour="black"),
                                                                                                                               axis.ticks = element_line(),
                                                                                                                               legend.title = element_text(face="italic"),
                                                                                                                               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p12



try<-data.frame(result4)
try = try[try$i == 11, ]
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))

p13<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zi_a_t[11], linetype = "dashed", color = "black")
p13<- p13 +  labs(title = expression(paste("Parameter ", zeta["11, H"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                    size = rel(1.2), hjust = 0.5),
                                                                                                                                
                                                                                                                                axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                                axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                                axis.title.x = element_text(vjust = -0.2),
                                                                                                                                axis.line = element_line(colour="black"),
                                                                                                                                axis.ticks = element_line(),
                                                                                                                                legend.title = element_text(face="italic"),
                                                                                                                                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p13




colnames(result5) <-c("X1","X2", "i", "it")
result5<-data.frame(result5)


try<-data.frame(result5)
try = try[try$i == 1, ]
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))


p14<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zi_b_t[1], linetype = "dashed", color = "black")
p14<- p14 +  labs(title = expression(paste("Parameter ", zeta["1, L"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                   size = rel(1.2), hjust = 0.5),
                                                                                                                               
                                                                                                                               axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                               axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                               axis.title.x = element_text(vjust = -0.2),
                                                                                                                               axis.line = element_line(colour="black"),
                                                                                                                               axis.ticks = element_line(),
                                                                                                                               legend.title = element_text(face="italic"),
                                                                                                                               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p14



try<-data.frame(result5)
try = try[try$i == 11, ]
try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))

p15<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  + geom_hline(yintercept  = zi_b_t[11], linetype = "dashed", color = "black")
p15<- p15 +  labs(title = expression(paste("Parameter ", zeta["11, L"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                                                                                                                    size = rel(1.2), hjust = 0.5),
                                                                                                                                
                                                                                                                                axis.title = element_text(face = "bold",size = rel(1)),
                                                                                                                                axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                                axis.title.x = element_text(vjust = -0.2),
                                                                                                                                axis.line = element_line(colour="black"),
                                                                                                                                axis.ticks = element_line(),
                                                                                                                                legend.title = element_text(face="italic"),
                                                                                                                                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p15



mc3<-(p12|p14)/(p13|p15)  + plot_annotation(tag_levels = 'A')
ggsave(mc3, filename = "Figures/Simulation/FigureD2.pdf", units = "cm", width = 16*2, height = 9*2 )



result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)




p12<-ggplot(result8) + geom_line(aes(x = ite , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = sa_mean), col = "black")  + geom_hline(yintercept  = var(zi_a_ref[,1]), linetype = "dashed", color = "black")
p12<- p12 +  labs(title = expression(paste("Parameter ", sigma[H]^2)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                                                                                                                 size = rel(1.2), hjust = 0.5),
                                                                                                                             axis.title = element_text(face = "italic",size = rel(1)),
                                                                                                                             axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                             axis.title.x = element_text(vjust = -0.2),
                                                                                                                             axis.line = element_line(colour="black"),
                                                                                                                             axis.ticks = element_line(),
                                                                                                                             legend.title = element_text(face="italic"),
                                                                                                                             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p12


p13<-ggplot(result8) + geom_line(aes(x = ite , y = X2), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = sb_mean), col = "black")  + geom_hline(yintercept  = var(zi_b_ref[,1]), linetype = "dashed", color = "black")
p13<- p13 +  labs(title = expression(paste("Parameter ", sigma[L]^2)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                                                                                                                 size = rel(1.2), hjust = 0.5),
                                                                                                                             axis.title = element_text(face = "italic",size = rel(1)),
                                                                                                                             axis.title.y = element_text(angle=90,vjust =2),
                                                                                                                             axis.title.x = element_text(vjust = -0.2),
                                                                                                                             axis.line = element_line(colour="black"),
                                                                                                                             axis.ticks = element_line(),
                                                                                                                             legend.title = element_text(face="italic"),
                                                                                                                             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
p13


p<-(p13|p12)

mc4<-(p13|p12) + plot_annotation(tag_levels = 'A')
ggsave(mc4, filename = "Figures/Simulation/FigureD3.pdf", units = "cm", width = 16*2, height = 9*2 )


###########################
#####TABLE D.1 ############

rm(list = ls())

#load the libraries

list.of.packages <- c("LaplacesDemon")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


library(LaplacesDemon)
load("Code/04-Simulation/RESULT_Sim_50kite.RData")



mat_res_1<- matrix(0, 6,63)
mat_res_2<- matrix(0, 6,63)
mat_res_3<- matrix(0, 6,63)


list_mat<- list(mat_res_1, mat_res_2, mat_res_3)
list_seq<- list(seq(1, 50000, 1), seq(30001, 50000, 1), seq(30001, 50000, 10))

result1<- result[[1]]
result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)

result4<- result[[3]]
result4<-data.frame(result4)
colnames(result4) <-c("X1","X2", "i", "it")

result5<- result[[4]]
result5<-data.frame(result5)
colnames(result5) <-c("X1","X2", "i", "it")

result4_wide  = result4 %>%  spread(i, X1)
result5_wide  = result5 %>%  spread(i, X1)
# 

series_tab<-cbind(result1[,1:20], result3[,c(2,3,1)], 
                  result5_wide[,3:ncol(result5_wide)],
                  result4_wide[,3:ncol(result4_wide)])

for( i in 1:3){
  
  list_seq_x = list_seq[[i]]
  
  for( k in 1:63){
    
    series<-  series_tab[,k]
    
    acf_res<-acf(series[list_seq_x], 100, plot = F)
    list_mat[[i]][1,k] <- round(acf_res$acf[1+1],3)
    list_mat[[i]][2,k] <-round(acf_res$acf[10+1],3)
    list_mat[[i]][3,k] <-round(acf_res$acf[30+1],3)
    list_mat[[i]][4,k] <- round(1-  sum(duplicated(series[list_seq_x]))/length(series[list_seq_x])    , 3)
    
    ess<- ESS(series[list_seq_x])
    list_mat[[i]][5,k]= round(ess/length(series[list_seq_x]),3)
    
    
    q<-Geweke.Diagnostic(series[list_seq_x])
    if(q > 0){list_mat[[i]][6,k]= round(pnorm(q, lower.tail = F),3 )}else{ list_mat[[i]][6,k] = round(pnorm(q, lower.tail = T),3 ) }
    
    print(k)
  }
  
  print(i)
}
print("Metrics Generated!")


mat_res_1<- list_mat[[1]]
mat_res_2<- list_mat[[2]]
mat_res_3<- list_mat[[3]]



m1<-cbind(rowMeans(mat_res_1[,1:20]),  mat_res_1[,21],  mat_res_1[,22], mat_res_1[,23], rowMeans(mat_res_1[,24:43]), rowMeans(mat_res_1[,44:63]))
m2<-cbind(rowMeans(mat_res_2[,1:20]),  mat_res_2[,21],  mat_res_2[,22], mat_res_2[,23], rowMeans(mat_res_2[,24:43]), rowMeans(mat_res_2[,44:63]))
m3<-cbind(rowMeans(mat_res_3[,1:20]),  mat_res_3[,21],  mat_res_3[,22], mat_res_3[,23], rowMeans(mat_res_3[,24:43]), rowMeans(mat_res_3[,44:63]))

M<- data.frame(rbind(m1, m2, m3))
colnames(M)<- c("alpha", "gamma0", "gamma1", "phi", "zl" , "zh")
rownames(M)<- paste0(c("acf1", "acf10", "acf30", "Acc", "ESS", "CD"), "_", rep(1:3, each = 6))

#set the acceptance rate of zi_ik to 0.25. Centering affects the acceptance counting.
M[4,c(5,6)] <- c(0.25,0.25)
M[10,c(5,6)] <- c(0.25, 0.25)
M[16,c(5,6)]<- c(0.25,0.25)


write.csv2(M, file = "Figures/Simulation/TableD.1.csv" , row.names = T)
