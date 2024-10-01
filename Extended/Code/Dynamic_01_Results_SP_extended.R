############## Gibbs Sampler #######################
############## Preliminaries ####################### 

rm(list = ls())

#load the libraries

library(mvtnorm)
library(Rcpp)
library(RcppDist)
library(RcppParallel)
library(RcppArmadillo)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

args <- commandArgs(TRUE)

Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}

odd<-function(x){(1/(1+exp(-1*x)))} #logistic function
odd_inv<-function(x){ (1-(1/(1+exp(-1*x))))  } #

log_sum_exp = function(x) {  b<- max(x); return(b + log(sum(exp(x-b))))  } #define log sum exp function
softmax <- function (x) {exp(x - log_sum_exp(x))} #define softmax function


load("Data/Dynamic/DataEnv_SP_all.RData")
sourceCpp("Code/Model/Extended_MS_LS_FE.cpp")

a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )


for(i in 1:length(unique(EL_x$i) )){
  
  name = unique(EL_x$i)[i]
  resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
  plot.ts(resa$x, main = name)
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


EL_x$ith_jth<-paste0(EL_x$ith,"_",EL_x$jth)
EL_x$i_j<- EL_x$ith_jth
EL_x$t<-as.numeric(as.factor(EL_x$t))
DBplane$t<-as.numeric(as.factor(DBplane$t))
EL_princ$t<-as.numeric(as.factor(EL_princ$t))
DBplane$ith<-as.numeric(as.factor(DBplane$i))
DBplane$leaning<- (DBplane$leaning - (min(DBplane$leaning)-0.01))/((max(DBplane$leaning)+0.01) - (min(DBplane$leaning)-0.01)   )


############## SETUP 
pages_names<- c(unique(EL$i)) # a vector of page names

#Main Dimensions
N<- length(unique(EL_x$i_j)) #number of edges
M<- length(unique(EL_x$i_j))/2 #number of unique edges
Time<- length(unique(EL_x$t)) #number of unique dates
Npages<- length(pages_names) #number of unique pages
K <- as.numeric(args[1]) #number of states
D<- as.numeric(args[2])

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

#Initialization of Xi

xi <- data.frame(matrix(rep(1/K,Time*K), nrow =  Time, ncol = K , byrow=TRUE))
xi<-t(apply(xi,1, FUN = Multi ))
# 
# xi<-data.frame(xi)
# xi$t<-1:nrow(xi)
# #xi<-cbind(xi,as.numeric(unique(EL_princ$t)))
# colnames(xi)<-c("state1", "state2", "t")

#Initialization of Pi 
P = diag(K) 

for(i in 1:K){
  for(j in 1:K){
    
    if(i == j){
      
      P[i,j] = P[i,j] -  0.10
      
    }else{
      
      P[i,j] =  0.10/(K-1)
      
    }
  }}

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

########

K <- as.numeric(args[1]) #number of states
D<- as.numeric(args[2])
lam_ad_z = matrix(1.2 ,Npages, K)

mu_mat_z = list()
Sigma_ad_z = list()
zi = list()

for(k in 1:K){
  mu_mat_z[[k]]  = matrix(0, Npages, D)
  
  id<-diag(rep(1,D))
  Sigma_ad_zk <- rep(list(id),Npages)
  
  Sigma_ad_z[[k]] <- Sigma_ad_zk
  
  zi[[k]] = matrix(0, Npages, D)
}

Iterations<- 45000

start_time <- Sys.time()  
result = MCMC(princ_w = EL_princ$w, #edge weights from the network triangle dataset
              princ_ones = EL_princ$ones,  #vector of ones from the network triangle dataset
              x_w = EL_x$w, #edge weights from the network off-diagonal dataset
              x_ones =  EL_x$ones, #vector of ones from the network off-diagonal dataset
              leaning = DBplane$leaning, #text analysis leaning from the dbplane dataset
              lam_ad_beta = lam_ad_beta, #adaptive RW-MH lambda - individual effects
              mu_mat_beta = mu_mat_beta,  #adaptive RW-MH mu - individual effects
              Sigma_ad_beta = Sigma_ad_beta,  #adaptive RW-MH sigma - individual effects
              lam_ad_z = lam_ad_z, 
              mu_mat_z = mu_mat_z,  
              Sigma_ad_z = Sigma_ad_z, 
              beta = rnorm(Npages, 0 , 0.01), # starting value - individual effect
              xi_state = as.matrix(xi[, 1:K]),
              zi_in = zi,  
              mu_beta =0, #prior mean - individual effect
              sigma_beta = 15,  #prior sd - individual effect
              mu_z = rep(0, D) ,  #prior mean - latent coordinate - state a
              sigma_z = rep(10, K),  #starting sd - latent coordinate - state a
              s_a = 0.01,  #prior shape parameter of sigma (gamma(s_a, s_b))
              s_b = 0.01,  #prior rate parameter of sigma (gamma(s_a, s_b))
              phi = 50, #starting value phi
              gamma_0 = 0,  #starting value -  gamma_0
              gamma_1 = 0,  #starting value - gamma_1
              a_phi = a_phi, #prior shape parameter of phi (gamma(a_phi, b_phi))
              b_phi = b_phi,  #prior rate parameter of phi (gamma(a_phi, b_phi))
              a_gamma_0 = 0,   #prior mean parameter of gamma_0
              b_gamma_0 = 15,  #prior sd parameter of gamma_0
              a_gamma_1 = 0,    #prior mean parameter of gamma_1
              b_gamma_1 = 15,    #prior sd parameter of gamma_1
              omega = rep(1/K, K), 
              P = P, #transition Matrix P
              N = Npages, # Number of nodes
              D = D,
              K = K,
              Time = Time , #Number of times
              x_i = EL_x$ith,  #index vector of nodes in the network off-diagonal dataset
              DBplane_i  = DBplane$i,  #index vector of nodes in the leaning dataset
              prop_sd_gamma = 0.01, #proposal RWMH for the gamma parameters
              prop_sd_phi = 10,  #proposal RWMH for the phi parameter
              acc_beta = 0.25 , # target acceptance Adaptive RW-MH -  individual effects
              acc_zeta = 0.25 ,  # target acceptance Adaptive RW-MH -  latent coordinates - state a
              pivot = which(pages_names == "ABC.es"),  # pivot node (news outlet) to solve identification
              sign= 1,  # sign of the news outlet -1 left or 1 right
              interp_eq = 1, # 0-1 option for the use of the interpretation equation
              ms_eq = 1, # 0-1 option for the use of the markov switching
              rg_eq = 0, # 0-1 option random graph
              Iterations = Iterations #number of iterations
)
end_time <- Sys.time()
end_time-start_time

interval<-35000:45000

ResList = Map('*',result$logLik,result$xi_ite)
ResList = lapply(ResList, rowSums)
ResList = lapply(ResList, sum)
ResLik = unlist(ResList)

Eta = matrix(0, Time, K)

for(k in 1:K){
  
  subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
  subz_k_ref = as.matrix(subz_k[[length(interval)]])
  
  zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(subz_k_ref,x )}) #apply procrustes trasnformation
  zi_k = apply(simplify2array(zi_res), 1:2, mean)
  Eta[, k] = Likelihood_k(  EL_princ$w,  DBplane$leaning,   beta = colMeans(result$beta_it[interval,]),  zi_k,  gamma_0 = mean(result$phi_gamma_ite[interval,2]),  gamma_1 = mean(result$phi_gamma_ite[interval,3]),  phi = mean(result$phi_gamma_ite[interval,1]),  N = Npages,  Time, 1);
  
}

xis = apply(simplify2array(result$xi_ite[interval]), 1:2, median)


likeBar = sum(Eta*xis)

DIC<- -2*mean(ResLik[interval]) +2*(likeBar - mean(ResLik[interval]) )
DIC

result$DIC<- DIC 

save(result, file = paste0("Extended/RESULT_SC_SPz_fe_extended_K_",as.numeric(args[1]),"_D_",as.numeric(args[2]),".Rdata"))