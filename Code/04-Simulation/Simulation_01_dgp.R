

list.of.packages <- c("data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")

library(data.table)

########CHANGE YOUR PATH ###########
#if you run this code individually
# setwd("~/Desktop/Repository/")
#####################################

#simulation
rm(list=ls())
set.seed(103)

Multi <- function(x){m<-rmultinom(1, size = 1, prob = x)
m<-t(m) 
return(m)
}

logistic<-function(x){(1/(1+exp(-1*x)))} #logistic function

Time<-100 # time
Npages<-20 #number of pages
N<-Npages*(Npages-1) # off diagonal elements
M<-N/2 # triangle

beta =  0 + rnorm(20,0, 2) 

zi_a_g1<-matrix(c(0, 0), nrow = Npages/2, ncol = 2,  byrow = T)
zi_a_g2<-matrix(c(3, 0), nrow = Npages/2, ncol = 2,  byrow = T)
zi_a<-rbind(zi_a_g1, zi_a_g2)
zi_a[,1]<-zi_a[,1]- mean(zi_a[,1])
zi_a[,1]<- ( 0.5*zi_a[,1] +rnorm(20,0, 0.15))


zi_b_g1<-matrix(c(1, 3), nrow = Npages/2, ncol = 2,  byrow = T) 
zi_b_g2<-matrix(c(2, 3), nrow = Npages/2, ncol = 2,  byrow = T)
zi_b<-rbind(zi_b_g1, zi_b_g2)
zi_b[,1]<-zi_b[,1]- mean(zi_b[,1])
zi_b[,1]<- (0.5*zi_b[,1] +rnorm(20,0, 0.15))

m_1 <-  mean( c(zi_a[,1], zi_b[,1] ))
m_2 <-  mean( c(zi_a[,2], zi_b[,2] ))


zi_a[,1]<-  zi_a[,1] -mean(zi_a[,1])
zi_a[,2] <- 0

zi_b[,1]<-  zi_b[,1] -mean(zi_b[,1])
zi_b[,2] <- 0

zi_a_ref<- zi_a
zi_b_ref<- zi_b

load("Data/Simulation/EL_sim.Rdata") #carico df vuoto

#empty data.frame
EL$t<-as.numeric(as.factor(EL$t))

#keep necessary dimensions
EL<-EL[EL$i %in% 1:Npages, ]
EL<-EL[EL$j %in% 1:Npages, ]
EL<-EL[EL$t %in% 1:Time, ] 

#create a starting number of common users
EL_start<-EL[EL$t == 1,]
EL_start$w<-1 # weight, y_ijt
EL_start$t<-"0"

# 
# EL_start$w[1:25] <- 100
# EL_start$w[(400-25+1):400] <- 100

#attach it to the previous df
EL<-rbind(EL_start, EL)
EL$t<-as.numeric(as.factor(EL$t))

#transition probabilities
p1<- c(0.95, 0.05)
p2<- c(0.05, 0.95)

P<-rbind(p1, p2)

phi <-  200 #beta precision
gamma1 <-  0.5 #beta location parameter
gamma0 <-   -0.1  #beta location parameter
tau <- 100 #pois intensity parameter

#states data frame
xi<-data.frame(t = 1:(Time+1), S1 = 0, S2= 0)
xi[1, 2:3]<- c(1,0)

#simulate transition probabilities
for(t in 2:(Time+1)){
  p<-as.matrix(xi[t-1, 2:3])%*%P
  xi[t, 2:3] <- Multi(p)
  print(t)
}

#compute the distance between pages
distance_a<- as.matrix(dist(cbind(zi_a[,1], rep(0, Npages)), method = "euclidean", diag = F, p = 2))
distance_b<- as.matrix(dist(cbind(zi_b[,1], rep(0, Npages)), method = "euclidean", diag = F, p = 2))

#rename the dataset
EL<-data.table(EL)

EL$t<-as.numeric(EL$t)
colnames(EL)<-c("ith", "jth", "th", "w", "wl1")

#create storing
M<-list()
Mat<-matrix( 0, nrow=Npages,  ncol=Npages)

Lam<-list()
Lam_mat<-matrix( 0, nrow=Npages,  ncol=Npages)

setkeyv(EL, c("th", "ith"))


for(t in 2:(Time+1)){
  S = xi$S1[t] #current state 
  
  for(i in 1:Npages){
    
    yt1<- unlist(EL[ith == i &  th == (t-1) , "w" ])
    
    lambda = exp( beta[i]  + beta   -1*(distance_a[,i]^2)*S  -1*(distance_b[,i]^2)*(1-S) ) 
    Lam_mat[i,]<- lambda
    aux <- rpois(Npages, lambda =  lambda) 
    EL[ith == i &  th == t , w := aux ] 
    Mat[i, ]<- aux
  }
  
  M[[t]]<-Mat
  Lam[[t]] <- Lam_mat
  
  print(t)}


# # 
EL$i_j<-paste0(EL$ith,"_",EL$jth)
EL$j_i<-paste0(EL$jth,"_",EL$ith)


for(i in unique(EL$j_i) ){
  EL[EL$i_j == i, "w"]<-EL[EL$j_i == i,"w"]
  
  print(i)
}

EL<-EL[, -7]

EL<-EL[order(EL$th), ]

setkeyv(EL,  c("th", "ith", "jth"))

for(i in unique(EL$ith)){
  for( j in unique(EL$jth)){
    
    w<- unlist(EL[ith == i & jth == j, "w" ])
    EL[ith == i & jth == j, wl1 :=   c(NA, w[1:Time] ) ] 
    print(paste0(i,"-", j))
  }}

colnames(EL)<- c("i","j","t" ,"w", "wl1", "i_j")
EL<-EL[EL$t != 1,]
EL$t<-as.numeric(as.factor(EL$t))
EL<-data.frame(EL)
EL<-EL[EL$i != EL$j,]
EL$ones<-1

EL<-EL[order(EL$t, EL$i, EL$j),]

EL_princ<- EL
EL_princ$i_t<-paste0(EL_princ$i,"_",EL_princ$t)
EL_princ$j_t<-paste0(EL_princ$j,"_",EL_princ$t)
EL_princ<-data.table(EL_princ)

setkeyv(EL_princ, c("i_j"))

A<-matrix(0, Npages,Npages)
ui<-which(upper.tri(A, diag = F), arr.ind=T)

upper_index<- paste0(ui[,1],"_",ui[,2])

EL_princ <- EL_princ[EL_princ$i_j %in% upper_index,    ]

EL_princ$ones<-1
EL_princ<-data.frame(EL_princ)

EL_princ<-EL_princ[, -c(8,9)]

EL_princ<-EL_princ[order(EL_princ$t, EL_princ$i, EL_princ$j),]


load("Data/Simulation/DBplane_sim.Rdata")

DBplane<-data.table(DBplane)

colnames(DBplane)[1]<-"ith"
colnames(DBplane)[4]<-"th"

DBplane$th<-as.numeric(as.factor(DBplane$th))
DBplane<-DBplane[DBplane$ith %in% 1:Npages]
DBplane<-DBplane[DBplane$th %in% 1:Time]

setkey(DBplane, th)

DBplane<-DBplane[order(DBplane$th,DBplane$ith )]


xi<-xi[-1, c(2,3,1)]
xi$t<-1:length(xi$t)

for(t in 1:Time ){
  
  S = xi$S1[t]
  
  lambda= tau*exp(S*zi_a[,2] + (1-S)*zi_b[,2])
  DBplane[th == t, ncomments:= unlist(lapply(lambda, FUN = function(x){rpois(1,x)}))] 
  
  
  beta_a <- phi*logistic(gamma0 + gamma1*(S*zi_a[,1]  + (1-S)*zi_b[,1] ) )
  beta_b <-phi*(1 - logistic(gamma0 + gamma1*(S*zi_a[,1] + (1-S)*zi_b[,1])  ) )
  beta_par<-data.frame(beta_a = beta_a, beta_b = beta_b)
  DBplane[th == t, leaning := apply(beta_par, 1, FUN = function(x){rbeta(1, shape1 = x[1], shape2 = x[2])})] 
  
  print(t)
}


colnames(DBplane)<-c("i","ncomments","leaning","t")
DBplane<-DBplane[order(DBplane$t,DBplane$i),]

#save(DBplane, file = "DBplane_sim_res.RData")

zi_a_t = zi_a[,1]
zi_b_t = zi_b[,1]
beta_t = beta  

remove_list<- ls()
remove_list<-remove_list[!remove_list %in% c("DBplane", "EL", "EL_princ","xi", "zi_a_ref", "zi_b_ref", "zi_a", "zi_b", "beta", "zi_a_t", "zi_b_t", "beta_t")]
rm(list= remove_list)



save.image("Data/Simulation/SimulationEnv_FE.RData")

print("1 - Synthetic Dataset Generated!")
