### Uncomment and run only
## if Results.Rdata files have been generate by Dynamic_01_Results_(Country)_extended.R
## Modify the script accordingly where written "Change...".

# library(mvtnorm)
# library(Rcpp)
# library(RcppDist)
# library(RcppParallel)
# library(RcppArmadillo)
# library(ggplot2)
# library(ggrepel)
# library(tidyr)
# library(scales)
# library(dplyr)
# library(data.table)
# library(patchwork)
# library(ggpubr)
# library(ggpmisc)
# library(lubridate)
# library(MCMCpack)
# library(gtable)
# library(cowplot)
# library(grid)
# 
# interval = 25000:45000
# 
# ########CHANGE YOUR PATH ###########
# setwd("~/Repository/")
# #####################################
# 
# load("Data/Dynamic/DataEnv_FR_all.RData")
# 
# ######## LOAD RESULTS  ############
# # load results.Rdata from Dynamic_01_Results_FR_extended.R
# load("")
# ###################################
# 
# sourceCpp("Code/Model/Extended_MS_LS_FE.cpp")
# sourceCpp("Code/Model/Extended_Predictive.cpp")
# 
# a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )
# 
# for(i in 1:length(unique(EL_x$i) )){
#   
#   name = unique(EL_x$i)[i]
#   resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
#   x<- rle(resa$x==0)
#   
#   if(length(x$lengths[x$values == TRUE])>0){
#     a[i,2]= max(x$lengths[x$values == TRUE])}
#   print(i)
# }
# 
# 
# thr <- c(a$name[a$max_czeros > 15])
# 
# 
# DBplane<-DBplane[! DBplane$fb_name %in% thr ,  ]
# DBplane$i<- as.numeric(factor(DBplane$fb_name))
# #
# EL_x<-EL_x[! EL_x$i %in% thr ,  ]
# EL_x<-EL_x[! EL_x$j %in% thr ,  ]
# #
# EL_princ<-EL_princ[! EL_princ$i %in% thr ,  ]
# EL_princ<-EL_princ[! EL_princ$j %in% thr ,  ]
# #
# EL<-EL[! EL$i %in% thr ,  ]
# EL<-EL[! EL$j %in% thr ,  ]
# 
# EL$ith<-as.numeric(factor(EL$i) , levels = unique(EL$i))
# EL$jth<-as.numeric(factor(EL$j) , levels = unique(EL$i))
# 
# EL_x$ith<-as.numeric(factor(EL_x$i) , levels = unique(EL_x$i))
# EL_x$jth<-as.numeric(factor(EL_x$j) , levels = unique(EL_x$i))
# 
# EL_princ$ith<-as.numeric(factor(EL_princ$i , levels = unique(EL_x$i)))
# EL_princ$jth<-as.numeric(factor(EL_princ$j, levels = unique(EL_x$i)))
# 
# pages_names_fr<- c(unique(EL_x$i)) # a vector of page names
# 
# #Main Dimensions
# N<- length(unique(EL_x$i_j)) #number of edges
# M<- length(unique(EL_x$i_j))/2 #number of unique edges
# Time<- length(unique(EL_x$t)) #number of unique dates
# Npages<- length(pages_names_fr) #number of unique pages
# K <- 2 #number of states
# D<-2
# 
# dates<- unique(as_date(DBplane$t))
# 
# 
# top_names<-c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart")
# 
# beta_fr<- result$beta_it
# agg_beta_fr<-colMeans(beta_fr[interval, ])
# 
# #################
# k = K
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# #################
# k = 1
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_1 = rbindlist(zi_res)
# zi_res_long_1 = data.frame(zi_res_long_1)
# zi_res_long_1$i = 1:Npages
# zi_res_long_1$it<- rep(1:length(zi_res), each = Npages)
# 
# agg_zeta_fr_a= aggregate(zi_res_long_1[,1:D], by = list(i= zi_res_long_1$i), FUN = mean)
# agg_zeta_fr_a$names = pages_names_fr
# colnames(agg_zeta_fr_a)[2:(D+1)]<-paste0("X",1:D)
# 
# #################
# k = 2
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# # subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_2 = rbindlist(zi_res)
# zi_res_long_2 = data.frame(zi_res_long_2)
# zi_res_long_2$i = 1:Npages
# 
# agg_zeta_fr_b= aggregate(zi_res_long_2[,1:D], by = list(i= zi_res_long_2$i), FUN = mean)
# agg_zeta_fr_b$names = pages_names_fr
# colnames(agg_zeta_fr_b)[2:(D+1)]<-paste0("X",1:D)
# 
# ############
# 
# if(K == 3){
#   k = 3
#   subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
#   subz_k_ref = as.matrix(subz_k[[length(interval)]])
#   
#   zi_res<-subz_k
#   zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
#   zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
#   
#   zi_res_long_3 = rbindlist(zi_res)
#   zi_res_long_3 = data.frame(zi_res_long_3)
#   zi_res_long_3$i = 1:Npages
#   
#   agg_zeta_fr_c= aggregate(zi_res_long_3[,1:D], by = list(i= zi_res_long_3$i), FUN = mean)
#   agg_zeta_fr_c$names = pages_names_fr
#   colnames(agg_zeta_fr_c)[2:(D+1)]<-paste0("X",1:D)
#   
# }
# 
# ##########
# 
# pew_fr_a<- agg_zeta_fr_a[agg_zeta_fr_a$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
# pew_fr_a$State<- "State H"
# pew_fr_a$score_a<- c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
# pew_fr_a$score_b<- c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3) - mean(c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3))
# 
# cor(pew_fr_a$X1, pew_fr_a$score_a)
# 
# pew_fr_b<- agg_zeta_fr_b[agg_zeta_fr_b$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
# pew_fr_b$State<- "State L"
# pew_fr_b$score_a<-  c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
# pew_fr_b$score_b<-  c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3) - mean(c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3))
# 
# pew_fr<-rbind(pew_fr_a, pew_fr_b)
# 
# 
# if(K == 3){
#   pew_fr_b$State<- "State M"
#   pew_fr_c<- agg_zeta_fr_c[agg_zeta_fr_c$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
#   pew_fr_c$State<- "State L"
#   pew_fr_c$score_a<-  c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
#   pew_fr_c$score_b<-  c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3) - mean(c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3))
#   pew_fr<-rbind(pew_fr_a, pew_fr_b, pew_fr_c)
# }
# 
# pew_fr$country<- "France"
# 
# db_fr_a<-data.frame(agg_zeta_fr_a[,2:(D+1)])
# colnames(db_fr_a)   <- paste0("agg_zeta_fr_", 1:D)             
# db_fr_a$agg_beta_fr= agg_beta_fr
# db_fr_a$names=pages_names_fr 
# db_fr_a$State<- "State H"
# 
# db_fr_b<-data.frame(agg_zeta_fr_b[,2:(D+1)])
# colnames(db_fr_b)   <- paste0("agg_zeta_fr_", 1:D)             
# db_fr_b$agg_beta_fr= agg_beta_fr
# db_fr_b$names=pages_names_fr 
# db_fr_b$State<- "State L"
# 
# if(K == 2){
#   db_fr<- rbind(db_fr_a,db_fr_b)
#   db_fr$State<-factor(db_fr$State,  levels = c("State L",  "State H"))
# }
# 
# if(K == 3){
#   db_fr_b$State<- "State M"  
#   
#   db_fr_c<-data.frame(agg_zeta_fr_c[,2:(D+1)])
#   colnames(db_fr_c)   <- paste0("agg_zeta_fr_", 1:D)             
#   db_fr_c$agg_beta_fr= agg_beta_fr
#   db_fr_c$names=pages_names_fr 
#   
#   db_fr_c$State<- "State L"
#   db_fr<- rbind(db_fr_a,db_fr_b, db_fr_c)
#   db_fr$State<-factor(db_fr$State,  levels = c("State L", "State M", "State H"))
# }
# 
# p2<-ggplot(db_fr, aes(x = agg_zeta_fr_1, y = agg_zeta_fr_2,  label = names)) +facet_wrap(~State)
# p2<- p2+ geom_point(col = "#ff420f", aes( size = exp(agg_beta_fr),))
# p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
# #p2<- p2 + geom_text_repel(data = db_fr, size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Latent Coordinate", title = "France", col = "white")
# #p2<- p2 + geom_text_repel(data = db_fr[db_fr$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "France", col = "white", size = "Individual Effect")
# p2<- p2 + geom_text_repel(data = db_fr,size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "France", col = "white", size = "Individual Effect")
# p2 <- p2 + theme(panel.grid = element_blank())  
# p2 <- p2 + theme_minimal() + theme(strip.placement = "outside",legend.position="none", plot.title = element_text(face = "italic", size = 15),
#                                    
#                                    
#                                    axis.title = element_text( face = "italic",size = rel(1)),
#                                    axis.title.y = element_text(size = 14, angle=90,vjust =2),
#                                    axis.title.x = element_blank(),
#                                    axis.text=element_text(size=12),
#                                    axis.line = element_blank(),
#                                    axis.ticks = element_line(),
#                                    panel.grid = element_blank(),
#                                    panel.border = element_rect(colour = "black", fill = NA),
#                                    strip.text.x = element_text(size = 14,face = "italic")  )
# 
# 
# 
# p2_fr<-p2
# 
# gp <- ggplotGrob(p2_fr)
# 
# 
# 
# 
# result3<- result$phi_gamma_ite
# 
# colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
# result3<-data.frame(result3)
# 
# 
# gat1<-gather(result3[result3$ite %in% interval,1:4], key = "parameter", value = "value")
# 
# result8<-data.frame(result$sigma_ite)
# result8$X1<- result8$X1^2
# result8$X2<- result8$X2^2
# if(K == 3){ result8$X3<- result8$X3^2}
# # 
# 
# result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
# result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
# if(K == 3){result8$sc_mean <- cumsum(result8$X3)/(1:length(result8$X3))}
# 
# result8$ite<- 1:length(result8$X2)
# 
# 
# gat2<-gather(result8[result8$ite %in% interval,1:K], key = "parameter", value = "value")
# 
# if(K == 2){gat2$parameter<-factor(gat2$parameter, levels=c("X2","X1"), labels = c("sigma[L]^2","sigma[H]^2"))
# }
# 
# if(K == 3){gat2$parameter<-factor(gat2$parameter, levels=c("X3", "X2","X1"), labels = c("sigma[L]^2","sigma[M]^2","sigma[H]^2"))
# }
# 
# 
# gat1$parameter<-factor(gat1$parameter)
# gat1$parameter<-factor(gat1$parameter, labels = c("gamma[0]", "gamma[1]", "phi", "tau"))
# 
# gat<-rbind(gat1, gat2)
# gat<-gat[!gat$parameter %in% c("alpha", "beta" , "tau"), ]
# 
# make_label <- function(value) {
#   x <- as.character(value)
#   bquote(italic(.(x))~subjects)
# }
# 
# plot_labeller <- function(variable, value) {
#   do.call(expression, lapply(levels(value), make_label))
# }
# 
# 
# gat<-gat[,1:2]
# 
# colnames(gat)[1]<-"level"
# 
# gat$type <- "posterior"
# gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[M]^2", "sigma[H]^2"))
# 
# gat_fr<- gat
# gat_fr$country <- "France"
# 
# sturges <- function(x){ pretty(range(x),
#                                n = nclass.Sturges(x),
#                                min.n = 1)}
# 
# z <- ggplot(gat,  aes(x=value ))
# z <-  z+ geom_histogram(data = gat[gat$level == "gamma[0]",], aes( y = after_stat(density)), binwidth=0.1, fill = "#ff420f", alpha = 0.6 , col ="black", breaks = sturges(gat[gat$level == "gamma[0]","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "gamma[1]",],  aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "gamma[1]","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "phi",],   aes( y = after_stat(density)),  fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "phi","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "sigma[L]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[L]^2","value"] ) ) 
# if(K == 3){ z <- z+ geom_histogram(data = gat[gat$level == "sigma[M]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[M]^2","value"] ) ) }
# z <- z+ geom_histogram(data = gat[gat$level == "sigma[H]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[H]^2","value"] ) ) 
# z<- z + geom_density(linetype = 1,  n = 10000, adjust = 2)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "phi") , fun = dgamma, args = list(shape=0.01, rate =0.01), geom = "area", fill ="black", alpha = 0.3, linetype = "dashed", colour = "black", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[0]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3, linetype = "dashed", colour = "black", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[1]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[L]^2"),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
# if(K == 3){z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[M]^2"),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)}
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[H]^2"),fun = dinvgamma, args = list(shape = 0.01, scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
# if(K == 3){z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[M]^2", "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed) }
# if(K == 2){z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2",  "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed) }
# z<- z + scale_shape_manual(labels =  parse_format())
# z <- z + labs(title = "France", 
#               x = "Value",  y = "")+  theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                 size = 22, hjust = 0.5),
#                                                               axis.text=element_text(size=8),
#                                                               strip.text.x = element_text(size = 18, face = "italic"),
#                                                               axis.title.x = element_blank(),
#                                                               axis.title.y = element_text(size = 22, angle=90,vjust =2, face = "italic"),
#                                                               axis.line = element_blank(),
#                                                               axis.ticks = element_line(),
#                                                               axis.ticks.x  = element_line(),
#                                                               panel.spacing.x = unit(4, "mm"),
#                                                               panel.grid = element_blank(),
#                                                               panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# w_fr<-z 
# 
# result7<- result$xi_ite
# result7<-lapply(result7, as.data.frame)
# result7<- data.frame(rbindlist(result7))
# result7$t <- 1:Time
# result7$ite<- rep(1:45000, each = Time)
# 
# colnames(result7)<-c(paste0("state",1:K), "t", "ite")
# result7<-data.frame(result7)
# 
# 
# xi_fr_sub<-result7[result7$ite %in% seq(35000, 45000, 10),]
# xi_fr_agg<-aggregate(xi_fr_sub[,1:K], by = list(t = xi_fr_sub$t), FUN = mean)
# xi_fr_agg$State<-  apply(xi_fr_agg[,2:(K+1)],1, which.max)
# # xi_it_agg$State<- factor(xi_it_agg$state1)
# # xi_it_agg$State<- factor(xi_it_agg$State, labels = c("L", "H"))
# xi_fr_agg$t <- dates
# 
# if(K == 2){
#   xi_fr_agg$State<-factor((K+1)-xi_fr_agg$State, labels = c("L", "H"))}
# if(K == 3){xi_fr_agg$State<-factor((K+1)-xi_fr_agg$State, labels = c("L", "M", "H"))}
# 
# q<- ggplot(xi_fr_agg)+labs(x ="time", y = "State", colour = "State")
# q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
# q<- q + geom_rect(data = xi_fr_agg , aes(xmin = t, xmax = lead(t),  ymax  = State),  ymin= "L", fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
# # q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
# #q<- q + geom_step(aes(x = t, y = (State), group = 1), color ="black", shape = 15 , size = 0.2)
# 
# # q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
# #                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# # 
# q<- q + labs(title = "France",
#              x = "Time", y = "State")+ theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                    size = 36, hjust = 0.5),
#                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                
#                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                axis.text=element_text(size=16),
#                                                                axis.line = element_blank(),
#                                                                axis.ticks = element_line(),
#                                                                panel.grid = element_blank(),
#                                                                legend.title = element_text(face="italic", size=16),
#                                                                legend.text = element_text(size=14),
#                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# q_fr<- q
# 
# 
# result6<- data.frame(result$P_ite)
# result6<- result6[seq(35000, 45000, 10),1:(K*K)]
# colnames(result6)<-paste("p",rep(1:K, each = length(1:K)), 1:K, sep = "")
# ps<- apply(result6,2, FUN = mean)
# ps
# ps_m<-matrix(ps, K,K, byrow = T)
# ps_m<-round(ps_m,3)
# 
# 
# # melted_cormat$Var1 <- relevel(melted_cormat$Var1, paste0(K))
# 
# if(K == 2){
#   
#   rownames(ps_m)<- c("H", "L")
#   colnames(ps_m)<- c("H", "L")
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L", "H"))
#   
#   }
# 
# if(K == 3){
#   
#   rownames(ps_m)<- c("H", "M", "L")
#   colnames(ps_m)<- c("H", "M", "L")
#   
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H","M", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L","M", "H"))
#   }
# 
# # Create a ggheatmap
# ggheatmap1 <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
#   geom_tile(color = "white")+ labs(x= "States", y ="States", title = "France" )+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   theme_minimal() + 
#   geom_text(aes(Var1, Var2, label = value), color = "black", size = 8) + theme(legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                                  size = 22, hjust = 0.5),
#                                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                                
#                                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                                axis.text=element_text(size=16),
#                                                                                axis.line = element_blank(),
#                                                                                axis.ticks = element_line(),
#                                                                                panel.grid = element_blank(),
#                                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# ggheatmap_fr <- ggheatmap1
# 
# EY = rep(0, Time)
# Var = rep(0, Time)
# DI = rep(0, Time)
# 
# agg<-aggregate(EL_x$w, by = list( i = EL_x$i, t = EL_x$t), FUN = sum) # strength distr through time
# 
# agg3<-aggregate(agg$x, by = list(agg$t), FUN = mean)
# agg4<-aggregate(agg$x, by = list(agg$t), FUN = var)
# 
# EY  = agg3$x
# Var = agg4$x
# DI = agg4$x/agg3$x
# 
# lkMetrics<- logPredictiveScoreMS_metrics(EL_princ$w,  result$beta_it, result$zi_ite,   result$xi_ite , Npages,  Time, K, D, 35*10^3, 45*10^3 , 10)
# 
# ResList_plot_db_long_a<-gather(data.frame(lkMetrics$EY), key = "Time", value = "Value")
# ResList_plot_db_long_a$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_a$Metric<- "Expected Strength"
# 
# ResList_plot_db_long_b<-gather(data.frame(lkMetrics$Var), key = "Time", value = "Value")
# ResList_plot_db_long_b$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_b$Metric<- "Variance Strength"
# 
# ResList_plot_db_long_c<-gather(data.frame(lkMetrics$DI), key = "Time", value = "Value")
# ResList_plot_db_long_c$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_c$Metric<- "Dispersion Index"
# 
# ResList_plot_db_long<-rbind(ResList_plot_db_long_a, ResList_plot_db_long_b, ResList_plot_db_long_c)
# ResList_plot_db_long$Metric<-factor(ResList_plot_db_long$Metric, levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ) )
# 
# q<-ggplot(ResList_plot_db_long, aes(x = Time, y = Value, group = Time), col ="red") + facet_wrap( vars(Metric), scales = "free")
# #p <- p +  geom_step(data = line_db, aes(x = date, y = value))
# q <- q +  geom_line(data = data.frame(Time = dates , Value = EY ,  Metric =   factor("Expected Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) ,col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = Var ,  Metric =   factor("Variance Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) , col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = DI ,  Metric =   factor("Dispersion Index", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value), col ="grey",  inherit.aes = F)
# q <- q +  geom_boxplot( fill = "red", col ="darkred", outlier.shape  = NA, lwd = 0.1)
# # p <- p +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
# #                         geom="crossbar", width=0.5)
# q<- q+ labs(y = "", title = "France") + theme_bw()
# q<- q + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                   size = 22, hjust = 0.5),
#               strip.text.x = element_text(size = 16,face = "italic"),
#               
#               axis.title = element_text( face = "italic",size = 16),
#               axis.title.y = element_text(size = 18, angle=90,vjust =2),
#               axis.title.x = element_text(size = 18, vjust = -0.2),
#               axis.text=element_text(size=8),
#               axis.line = element_blank(),
#               axis.ticks = element_line(),
#               panel.grid = element_blank(),
#               legend.title = element_text(face="italic", size=16),
#               legend.text = element_text(size=14),
#               panel.border = element_rect(colour = "black", fill = NA))
# gv_fr<-q
# 
# ##########################
# setwd("~/Repository/")
# #####################################
# 
# load("Data/Dynamic/DataEnv_DE_all.RData")

# ######## LOAD RESULTS  ############
# # load results.Rdata from Dynamic_01_Results_DE_extended.R
# load("")
# ###################################
# 
# sourceCpp("Code/Model/Extended_MS_LS_FE.cpp")
# sourceCpp("Code/Model/Extended_Predictive.cpp")
# 
# 
# a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )
# 
# 
# for(i in 1:length(unique(EL_x$i) )){
#   
#   name = unique(EL_x$i)[i]
#   resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
#   x<- rle(resa$x==0)
#   
#   if(length(x$lengths[x$values == TRUE])>0){
#     a[i,2]= max(x$lengths[x$values == TRUE])}
#   print(i)
# }
# 
# 
# thr <- c(a$name[a$max_czeros > 15])
# 
# 
# DBplane<-DBplane[! DBplane$fb_name %in% thr ,  ]
# DBplane$i<- as.numeric(factor(DBplane$fb_name))
# #
# EL_x<-EL_x[! EL_x$i %in% thr ,  ]
# EL_x<-EL_x[! EL_x$j %in% thr ,  ]
# #
# EL_princ<-EL_princ[! EL_princ$i %in% thr ,  ]
# EL_princ<-EL_princ[! EL_princ$j %in% thr ,  ]
# #
# EL<-EL[! EL$i %in% thr ,  ]
# EL<-EL[! EL$j %in% thr ,  ]
# 
# EL$ith<-as.numeric(factor(EL$i) , levels = unique(EL$i))
# EL$jth<-as.numeric(factor(EL$j) , levels = unique(EL$i))
# 
# EL_x$ith<-as.numeric(factor(EL_x$i) , levels = unique(EL_x$i))
# EL_x$jth<-as.numeric(factor(EL_x$j) , levels = unique(EL_x$i))
# 
# EL_princ$ith<-as.numeric(factor(EL_princ$i , levels = unique(EL_x$i)))
# EL_princ$jth<-as.numeric(factor(EL_princ$j, levels = unique(EL_x$i)))
# 
# pages_names_de<- c(unique(EL_x$i)) # a vector of page names
# pages_names_de[10]<- "FAZ"
# 
# #Main Dimensions
# N<- length(unique(EL_x$i_j)) #number of edges
# M<- length(unique(EL_x$i_j))/2 #number of unique edges
# Time<- length(unique(EL_x$t)) #number of unique dates
# Npages<- length(pages_names_de) #number of unique pages
# K <- 2 #number of states
# D<-2
# 
# dates<- unique(as_date(EL_princ$tindex))
# 
# 
# top_names<-  c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" )
# 
# beta_de<- result$beta_it
# agg_beta_de<-colMeans(beta_de[interval, ])
# 
# #################
# k = 1
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# #################
# k = 1
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# # subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_1 = rbindlist(zi_res)
# zi_res_long_1 = data.frame(zi_res_long_1)
# zi_res_long_1$i = 1:Npages
# zi_res_long_1$it<- rep(1:length(zi_res), each = Npages)
# 
# agg_zeta_de_a= aggregate(zi_res_long_1[,1:D], by = list(i= zi_res_long_1$i), FUN = mean)
# agg_zeta_de_a$names = pages_names_de
# colnames(agg_zeta_de_a)[2:(D+1)]<-paste0("X",1:D)
# 
# #################
# k = 2
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# # subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_2 = rbindlist(zi_res)
# zi_res_long_2 = data.frame(zi_res_long_2)
# zi_res_long_2$i = 1:Npages
# 
# agg_zeta_de_b= aggregate(zi_res_long_2[,1:D], by = list(i= zi_res_long_2$i), FUN = mean)
# agg_zeta_de_b$names = pages_names_de
# colnames(agg_zeta_de_b)[2:(D+1)]<-paste0("X",1:D)
# 
# ############
# 
# if(K == 3){
#   k = 3
#   subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
#   subz_k_ref = as.matrix(subz_k[[length(interval)]])
#   
#   zi_res<-subz_k
#   zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
#   zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
#   
#   zi_res_long_3 = rbindlist(zi_res)
#   zi_res_long_3 = data.frame(zi_res_long_3)
#   zi_res_long_3$i = 1:Npages
#   
#   agg_zeta_de_c= aggregate(zi_res_long_3[,1:D], by = list(i= zi_res_long_3$i), FUN = mean)
#   agg_zeta_de_c$names = pages_names_de
#   colnames(agg_zeta_de_c)[2:(D+1)]<-paste0("X",1:D)
#   
#   agg_zeta_de_c$X1  = -1*agg_zeta_de_c$X1
#  
# }
# 
# ##############
# 
# pew_de_a<- agg_zeta_de_a[agg_zeta_de_a$names %in%c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
# pew_de_a$State<- "State H"
# pew_de_a$score_a<- c(3.6, 3.2, 3.2, 2.7, 3) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
# pew_de_a$score_b<- c(3.1, 2.9, 3, 2.8, 2.8) - mean(c(3.1, 2.9, 3, 2.8, 2.8))
# cor(pew_de_a$X1, pew_de_a$score_a)
# 
# pew_de_b<- agg_zeta_de_b[agg_zeta_de_b$names %in%  c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
# pew_de_b$State<- "State L"
# pew_de_b$score_a<- c(3.6, 3.2, 3.2, 2.7, 3) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
# pew_de_b$score_b<- c(3.1, 2.9, 3, 2.8, 2.8) - mean(c(3.1, 2.9, 3, 2.8, 2.8))
# 
# if(K == 3){
#   pew_de_b$State<- "State M"
#   pew_de_c<- agg_zeta_de_c[agg_zeta_de_c$names %in%  c("Bild", "FAZ", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
#   pew_de_c$State<- "State L"
#   pew_de_c$score_a<- c(3.6, 3.2, 3.2, 2.7, 3) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
#   pew_de_c$score_b<-  c(3.1, 2.9, 3, 2.8, 2.8) - mean(c(3.1, 2.9, 3, 2.8, 2.8))
#   pew_de<-rbind(pew_de_a, pew_de_b, pew_de_c)
# }
# 
# pew_de<-rbind(pew_de_a, pew_de_b)
# 
# pew_de$country<- "Germany"
# 
# db_de_a<-data.frame(agg_zeta_de_a[,2:(D+1)])
# colnames(db_de_a)   <- paste0("agg_zeta_de_", 1:D)             
# db_de_a$agg_beta_de= agg_beta_de
# db_de_a$names=pages_names_de 
# db_de_a$State<- "State H"
# 
# 
# db_de_b<-data.frame(agg_zeta_de_b[,2:(D+1)])
# colnames(db_de_b)   <- paste0("agg_zeta_de_", 1:D)             
# db_de_b$agg_beta_de= agg_beta_de
# db_de_b$names=pages_names_de
# db_de_b$State<- "State L"
# 
# db_de<- rbind(db_de_a,db_de_b)
# db_de$State<-factor(db_de$State,  levels = c("State L", "State H"))
# 
# 
# if(K == 3){
#   db_de_b$State<- "State M"
#   
#   db_de_c<-data.frame(agg_zeta_de_c[,2:(D+1)])
#   colnames(db_de_c)   <- paste0("agg_zeta_de_", 1:D)             
#   db_de_c$agg_beta_de= agg_beta_de
#   db_de_c$names=pages_names_de
#   
#   db_de_c$State<- "State L"
#   db_de<- rbind(db_de_a,db_de_b, db_de_c)
#   db_de$State<-factor(db_de$State,  levels = c("State L", "State M", "State H"))
# }
# 
# 
# p2<-ggplot(db_de, aes(x = agg_zeta_de_1, y = agg_zeta_de_2,  label = names)) +facet_wrap(~State)
# p2<- p2+ geom_point(col = "#ff420f", aes( size = exp(agg_beta_de),))
# p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
# #p2<- p2 + geom_text_repel(data = db_de, size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Latent Coordinate", title = "Germany", col = "white")
# #p2<- p2 + geom_text_repel(data = db_de[db_de$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "Germany", col = "white", size = "Individual Effect")
# p2<- p2 + geom_text_repel(data = db_de,size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "Germany", col = "white", size = "Individual Effect")
# p2 <- p2 + theme(panel.grid = element_blank())  
# p2 <- p2 + theme_minimal() + theme(strip.placement = "outside",legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                  size = 14),
#                                    strip.text.x = element_blank(),
#                                    
#                                    axis.title = element_text( face = "italic",size = rel(1)),
#                                    axis.title.y = element_text(size = 14, angle=90,vjust =2),
#                                    axis.title.x = element_blank(),
#                                    axis.text=element_text(size=12),
#                                    axis.line = element_blank(),
#                                    axis.ticks = element_line(),
#                                    panel.grid = element_blank(),
#                                    panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# 
# 
# p2_de<-p2
# 
# 
# result3<- result$phi_gamma_ite
# 
# colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
# result3<-data.frame(result3)
# 
# 
# gat1<-gather(result3[result3$ite %in% interval,1:4], key = "parameter", value = "value")
# 
# result8<-data.frame(result$sigma_ite)
# result8$X1<- result8$X1^2
# result8$X2<- result8$X2^2
# if(K == 3){result8$X3<- result8$X3^2}
# 
# result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
# result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
# if(K == 3){result8$sc_mean <- cumsum(result8$X3)/(1:length(result8$X3))}
# 
# result8$ite<- 1:length(result8$X2)
# 
# 
# gat2<-gather(result8[result8$ite %in% interval,1:K], key = "parameter", value = "value")
# if(K == 3){gat2$parameter<-factor(gat2$parameter, levels=c("X3", "X2","X1"), labels = c("sigma[L]^2","sigma[M]^2","sigma[H]^2"))}
# if(K == 2){gat2$parameter<-factor(gat2$parameter, levels=c("X2","X1"), labels = c("sigma[L]^2","sigma[H]^2"))}
# 
# 
# gat1$parameter<-factor(gat1$parameter)
# gat1$parameter<-factor(gat1$parameter, labels = c("gamma[0]", "gamma[1]", "phi", "tau"))
# 
# gat<-rbind(gat1, gat2)
# gat<-gat[!gat$parameter %in% c("alpha", "beta" , "tau"), ]
# 
# 
# 
# make_label <- function(value) {
#   x <- as.character(value)
#   bquote(italic(.(x))~subjects)
# }
# 
# plot_labeller <- function(variable, value) {
#   do.call(expression, lapply(levels(value), make_label))
# }
# 
# 
# gat<-gat[,1:2]
# 
# colnames(gat)[1]<-"level"
# 
# gat$type <- "posterior"
# if(K == 2){gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2",  "sigma[H]^2"))}
# if(K == 3){gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[M]^2",  "sigma[H]^2"))}
# 
# gat_de<- gat
# gat_de$country <- "Germany"
# 
# sturges <- function(x){ pretty(range(x),
#                                n = nclass.Sturges(x),
#                                min.n = 1)}
# 
# z <- ggplot(gat,  aes(x=value ))
# z <-  z+ geom_histogram(data = gat[gat$level == "gamma[0]",], aes( y = after_stat(density)), binwidth=0.2, fill = "#ff420f", alpha = 0.6 , col ="black", breaks = sturges(gat[gat$level == "gamma[0]","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "gamma[1]",],  aes( y = after_stat(density)), binwidth=0.02, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "gamma[1]","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "phi",],   aes( y = after_stat(density)), binwidth=0.02, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "phi","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "sigma[L]^2",],   aes( y = after_stat(density)), binwidth=0.02, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[L]^2","value"] ) ) 
# if(K == 3){ z <- z+ geom_histogram(data = gat[gat$level == "sigma[M]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[M]^2","value"] ) ) }
# z <- z+ geom_histogram(data = gat[gat$level == "sigma[H]^2",],   aes( y = after_stat(density)), binwidth=0.02, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[H]^2","value"] ) ) 
# z<- z + geom_density(linetype = 1,  n = 10000, adjust = 2)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "phi") , fun = dgamma, args = list(shape=0.01, rate =0.01), geom = "area", fill ="black", alpha = 0.3, linetype = "dashed", colour = "black", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[0]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3, linetype = "dashed", colour = "black", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[1]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[L]^2"),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
# if(K == 3){ z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[M]^2"),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)}
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[H]^2"),fun = dinvgamma, args = list(shape = 0.01, scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "black", inherit.aes = F)
# if(K == 3){ z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[M]^2", "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed) }
# if(K == 2){z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed) }
# z<- z + scale_shape_manual(labels =  parse_format())
# z <- z + labs(title = "Germany", 
#               x = "Value",  y = "")+  theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                 size = 22, hjust = 0.5),
#                                                               axis.text=element_text(size=8),
#                                                               strip.text.x = element_text(size = 18, face = "italic"),
#                                                               axis.title.x = element_blank(),
#                                                               axis.title.y = element_text(size = 22, angle=90,vjust =2, face = "italic"),
#                                                               axis.line = element_blank(),
#                                                               axis.ticks = element_line(),
#                                                               axis.ticks.x  = element_line(),
#                                                               panel.spacing.x = unit(4, "mm"),
#                                                               panel.grid = element_blank(),
#                                                               panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# w_de<-z 
# 
# 
# result7<- result$xi_ite
# result7<-lapply(result7, as.data.frame)
# result7<- data.frame(rbindlist(result7))
# result7$t <- 1:Time
# result7$ite<- rep(1:45000, each = Time)
# 
# colnames(result7)<-c(paste0("state",1:K), "t", "ite")
# result7<-data.frame(result7)
# 
# xi_de_sub<-result7[result7$ite %in% seq(35000, 45000, 10),]
# xi_de_agg<-aggregate(xi_de_sub[,1:K], by = list(t = xi_de_sub$t), FUN = mean)
# xi_de_agg$State<-  apply(xi_de_agg[,2:(K+1)],1, which.max)
# # xi_it_agg$State<- factor(xi_it_agg$state1)
# # xi_it_agg$State<- factor(xi_it_agg$State, labels = c("L", "H"))
# xi_de_agg$t <- dates
# if(K == 2){xi_de_agg$State<-factor((K+1)-xi_de_agg$State, labels = c("L", "H"))}
# if(K == 3){xi_de_agg$State<-factor((K+1)-xi_de_agg$State, labels = c("L","M", "H"))}
# 
# q<- ggplot(xi_de_agg)+labs(x ="time", y = "State", colour = "State")
# q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
# q<- q + geom_rect(data = xi_de_agg , aes(xmin = t, xmax = lead(t),  ymax  = State),  ymin= "L", fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
# # q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
# #q<- q + geom_step(aes(x = t, y = (State), group = 1), color ="black", shape = 15 , size = 0.2)
# # q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
# #                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# # 
# q<- q + labs(title = "Germany",
#              x = "Time", y = "State")+ theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                    size = 36, hjust = 0.5),
#                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                
#                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                axis.text=element_text(size=16),
#                                                                axis.line = element_blank(),
#                                                                axis.ticks = element_line(),
#                                                                panel.grid = element_blank(),
#                                                                legend.title = element_text(face="italic", size=16),
#                                                                legend.text = element_text(size=14),
#                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# q_de<- q
# 
# 
# 
# 
# result6<- data.frame(result$P_ite)
# result6<- result6[seq(35000, 45000, 10),1:(K*K)]
# colnames(result6)<-paste("p",rep(1:K, each = length(1:K)), 1:K, sep = "")
# ps<- apply(result6,2, FUN = mean)
# ps
# ps_m<-matrix(ps, K,K, byrow = T)
# ps_m<-round(ps_m,3)
# 
# 
# if(K == 2){
#   
#   rownames(ps_m)<- c("H", "L")
#   colnames(ps_m)<- c("H", "L")
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L", "H"))
#   
# }
# 
# if(K == 3){
#   
#   rownames(ps_m)<- c("H", "M", "L")
#   colnames(ps_m)<- c("H", "M", "L")
#   
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H","M", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L","M", "H"))
# }
# 
# # Create a ggheatmap
# ggheatmap1 <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
#   geom_tile(color = "white")+ labs(x= "States", y ="States", title = "Germany" )+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   theme_minimal() + 
#   geom_text(aes(Var1, Var2, label = value), color = "black", size = 8) + theme(legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                                  size = 22, hjust = 0.5),
#                                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                                
#                                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                                axis.text=element_text(size=16),
#                                                                                axis.line = element_blank(),
#                                                                                axis.ticks = element_line(),
#                                                                                panel.grid = element_blank(),
#                                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# ggheatmap_de <- ggheatmap1
# 
# EY = rep(0, Time)
# Var = rep(0, Time)
# DI = rep(0, Time)
# 
# agg<-aggregate(EL_x$w, by = list( i = EL_x$i, t = EL_x$t), FUN = sum) # strength distr through time
# 
# agg3<-aggregate(agg$x, by = list(agg$t), FUN = mean)
# agg4<-aggregate(agg$x, by = list(agg$t), FUN = var)
# 
# EY  = agg3$x
# Var = agg4$x
# DI = agg4$x/agg3$x
# 
# lkMetrics<- logPredictiveScoreMS_metrics(EL_princ$w,  result$beta_it, result$zi_ite,   result$xi_ite , Npages,  Time, K, D, 35*10^3, 45*10^3 , 10)
# 
# ResList_plot_db_long_a<-gather(data.frame(lkMetrics$EY), key = "Time", value = "Value")
# ResList_plot_db_long_a$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_a$Metric<- "Expected Strength"
# 
# ResList_plot_db_long_b<-gather(data.frame(lkMetrics$Var), key = "Time", value = "Value")
# ResList_plot_db_long_b$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_b$Metric<- "Variance Strength"
# 
# ResList_plot_db_long_c<-gather(data.frame(lkMetrics$DI), key = "Time", value = "Value")
# ResList_plot_db_long_c$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_c$Metric<- "Dispersion Index"
# 
# ResList_plot_db_long<-rbind(ResList_plot_db_long_a, ResList_plot_db_long_b, ResList_plot_db_long_c)
# ResList_plot_db_long$Metric<-factor(ResList_plot_db_long$Metric, levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ) )
# 
# q<-ggplot(ResList_plot_db_long, aes(x = Time, y = Value, group = Time), col ="red") + facet_wrap( vars(Metric), scales = "free")
# #p <- p +  geom_step(data = line_db, aes(x = date, y = value))
# q <- q +  geom_line(data = data.frame(Time = dates , Value = EY ,  Metric =   factor("Expected Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) ,col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = Var ,  Metric =   factor("Variance Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) , col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = DI ,  Metric =   factor("Dispersion Index", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value), col ="grey",  inherit.aes = F)
# q <- q +  geom_boxplot( fill = "red", col ="darkred", outlier.shape  = NA, lwd = 0.1)
# # p <- p +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
# #                         geom="crossbar", width=0.5)
# q<- q+ labs(y = "", title = "Germany") + theme_bw()
# q<- q + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                   size = 22, hjust = 0.5),
#               strip.text.x = element_text(size = 16,face = "italic"),
#               
#               axis.title = element_text( face = "italic",size = 16),
#               axis.title.y = element_text(size = 18, angle=90,vjust =2),
#               axis.title.x = element_text(size = 18, vjust = -0.2),
#               axis.text=element_text(size=8),
#               axis.line = element_blank(),
#               axis.ticks = element_line(),
#               panel.grid = element_blank(),
#               legend.title = element_text(face="italic", size=16),
#               legend.text = element_text(size=14),
#               panel.border = element_rect(colour = "black", fill = NA))
# gv_de<-q
# 
# ##########################
# 
# interval = 25000:45000
# 
# 
# list.of.packages <- c( "ggplot2", "ggrepel", "tidyr", "dplyr", "data.table", "patchwork", "ggpubr")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages, type = "binary")
# 
# 
# 
# interval = 25000:45000
# 
# ########CHANGE YOUR PATH ###########
# setwd("~/Repository/")
# #####################################
# 
# ######## LOAD RESULTS  ############
# # load results.Rdata from Dynamic_01_Results_IT_extended.R
# load("")
# ###################################
# 
# sourceCpp("Code/Model/Extended_MS_LS_FE.cpp")
# sourceCpp("Code/Model/Extended_Predictive.cpp")
# 
# 
# a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )
# 
# 
# for(i in 1:length(unique(EL_x$i) )){
#   
#   name = unique(EL_x$i)[i]
#   resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
#   x<- rle(resa$x==0)
#   
#   if(length(x$lengths[x$values == TRUE])>0){
#     a[i,2]= max(x$lengths[x$values == TRUE])}
#   print(i)
# }
# 
# 
# thr <- c(a$name[a$max_czeros > 15])
# 
# 
# DBplane<-DBplane[! DBplane$fb_name %in% thr ,  ]
# DBplane$i<- as.numeric(factor(DBplane$fb_name))
# #
# EL_x<-EL_x[! EL_x$i %in% thr ,  ]
# EL_x<-EL_x[! EL_x$j %in% thr ,  ]
# #
# EL_princ<-EL_princ[! EL_princ$i %in% thr ,  ]
# EL_princ<-EL_princ[! EL_princ$j %in% thr ,  ]
# #
# EL<-EL[! EL$i %in% thr ,  ]
# EL<-EL[! EL$j %in% thr ,  ]
# 
# EL$ith<-as.numeric(factor(EL$i) , levels = unique(EL$i))
# EL$jth<-as.numeric(factor(EL$j) , levels = unique(EL$i))
# 
# EL_x$ith<-as.numeric(factor(EL_x$i) , levels = unique(EL_x$i))
# EL_x$jth<-as.numeric(factor(EL_x$j) , levels = unique(EL_x$i))
# 
# EL_princ$ith<-as.numeric(factor(EL_princ$i , levels = unique(EL_x$i)))
# EL_princ$jth<-as.numeric(factor(EL_princ$j, levels = unique(EL_x$i)))
# 
# pages_names_it<- c(unique(EL_x$i)) # a vector of page names
# 
# #Main Dimensions
# N<- length(unique(EL_x$i_j)) #number of edges
# M<- length(unique(EL_x$i_j))/2 #number of unique edges
# Time<- length(unique(EL_x$t)) #number of unique dates
# Npages<- length(pages_names_it) #number of unique pages
# K <- 2 #number of states
# D<-2
# 
# dates<- unique(as_date(DBplane$t))
# 
# top_names<-c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" )
# 
# beta_it<- result$beta_it
# agg_beta_it<-colMeans(beta_it[interval, ])
# 
# #################
# k = K
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# #################
# k = 1
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# # subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_1 = rbindlist(zi_res)
# zi_res_long_1 = data.frame(zi_res_long_1)
# zi_res_long_1$i = 1:Npages
# zi_res_long_1$it<- rep(1:length(zi_res), each = Npages)
# 
# agg_zeta_it_a= aggregate(zi_res_long_1[,1:D], by = list(i= zi_res_long_1$i), FUN = mean)
# agg_zeta_it_a$names = pages_names_it
# colnames(agg_zeta_it_a)[2:(D+1)]<-paste0("X",1:D)
# 
# #################
# k = 2
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_2 = rbindlist(zi_res)
# zi_res_long_2 = data.frame(zi_res_long_2)
# zi_res_long_2$i = 1:Npages
# 
# agg_zeta_it_b= aggregate(zi_res_long_2[,1:D], by = list(i= zi_res_long_2$i), FUN = mean)
# agg_zeta_it_b$names = pages_names_it
# colnames(agg_zeta_it_b)[2:(D+1)]<-paste0("X",1:D)
# 
# ############
# 
# if(K == 3){
#   k = 3
#   subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
#   subz_k_ref = as.matrix(subz_k[[length(interval)]])
#   
#   zi_res<-subz_k
#   zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
#   zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
#   
#   zi_res_long_3 = rbindlist(zi_res)
#   zi_res_long_3 = data.frame(zi_res_long_3)
#   zi_res_long_3$i = 1:Npages
# 
#   agg_zeta_it_c= aggregate(zi_res_long_3[,1:D], by = list(i= zi_res_long_3$i), FUN = mean)
#   agg_zeta_it_c$names = pages_names_it
#   colnames(agg_zeta_it_c)[2:(D+1)]<-paste0("X",1:D)
#   
# }
# ##############
# 
# pew_it_a<- agg_zeta_it_a[agg_zeta_it_a$names %in%c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ), ]
# pew_it_a$State<- "State H"
# pew_it_a$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5) - mean( c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5))
# pew_it_a$score_b<-  c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7) - mean( c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)  )
# 
# cor(pew_it_a$X1, pew_it_a$score_a)
# 
# pew_it_b<- agg_zeta_it_b[agg_zeta_it_b$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ), ]
# pew_it_b$State<- "State L"
# pew_it_b$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5) - mean( c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5))
# pew_it_b$score_b<-  c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7) - mean( c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)  )
# pew_it<-rbind(pew_it_a, pew_it_b)
# 
# if(K == 3){
#   pew_it_b$State<- "State M"
#   pew_it_c<- agg_zeta_it_c[agg_zeta_it_c$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ), ]
#   pew_it_c$State<- "State L"
#   pew_it_c$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5) - mean( c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5))
#   pew_it_c$score_b<-  c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7) - mean( c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)  )
# }
# 
# if(K == 3){pew_it<-rbind(pew_it_a, pew_it_b, pew_it_c)}
# 
# pew_it$country<- "Italy"
# 
# db_it_a<-data.frame(agg_zeta_it_a[,2:(D+1)])
# colnames(db_it_a)   <- paste0("agg_zeta_it_", 1:D)             
# db_it_a$agg_beta_it= agg_beta_it
# db_it_a$names=pages_names_it 
# 
# db_it_a$State<- "State H"
# 
# db_it_b<-data.frame(agg_zeta_it_b[,2:(D+1)])
# colnames(db_it_b)   <- paste0("agg_zeta_it_", 1:D)             
# db_it_b$agg_beta_it= agg_beta_it
# db_it_b$names=pages_names_it 
# db_it_b$State<- "State L"
# 
# db_it<- rbind(db_it_a,db_it_b)
# db_it$State<-factor(db_it$State,  levels = c("State L", "State H"))
# 
# if(K == 3){
#   db_it_b$State<- "State M"
#   
#   db_it_c<-data.frame(agg_zeta_it_c[,2:(D+1)])
#   colnames(db_it_c)   <- paste0("agg_zeta_it_", 1:D)             
#   db_it_c$agg_beta_it= agg_beta_it
#   db_it_c$names=pages_names_it 
#   
#   db_it_c$State<- "State L"
#   
#   db_it<- rbind(db_it_a,db_it_b, db_it_c)
#   db_it$State<-factor(db_it$State,  levels = c("State L", "State M", "State H"))
#   
# }
# 
# 
# 
# p2<-ggplot(db_it, aes(x = agg_zeta_it_1, y = agg_zeta_it_2,  label = names)) +facet_wrap(~State)
# p2<- p2+ geom_point(col = "#ff420f", aes( size = exp(agg_beta_it),))
# p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
# #p2<- p2 + geom_text_repel(data = db_it, size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Latent Coordinate", title = "Italy", col = "white")
# #p2<- p2 + geom_text_repel(data = db_it[db_it$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "Italy", col = "white", size = "Individual Effect")
# p2<- p2 + geom_text_repel(data = db_it,size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "Italy", col = "white", size = "Individual Effect")
# p2 <- p2 + theme(panel.grid = element_blank())  
# p2 <- p2 + theme_minimal() + theme(strip.placement = "outside",legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                  size = 14),
#                                    strip.text.x = element_blank(),
#                                    
#                                    axis.title = element_text( face = "italic",size = rel(1)),
#                                    axis.title.y = element_text(size = 14, angle=90,vjust =2),
#                                    axis.title.x = element_blank(),
#                                    axis.text=element_text(size=12),
#                                    axis.line = element_blank(),
#                                    axis.ticks = element_line(),
#                                    panel.grid = element_blank(),
#                                    panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# 
# 
# p2_it<-p2
# 
# 
# result3<- result$phi_gamma_ite
# 
# colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
# result3<-data.frame(result3)
# 
# 
# gat1<-gather(result3[result3$ite %in% interval,1:4], key = "parameter", value = "value")
# 
# result8<-data.frame(result$sigma_ite)
# result8$X1<- result8$X1^2
# result8$X2<- result8$X2^2
# if(K == 3){ result8$X3<- result8$X3^2}
# 
# result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
# result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
# if(K == 3){ result8$sc_mean <- cumsum(result8$X3)/(1:length(result8$X3))}
# 
# result8$ite<- 1:length(result8$X2)
# 
# 
# gat2<-gather(result8[result8$ite %in% interval,1:K], key = "parameter", value = "value")
# if(K == 2){  gat2$parameter<-factor(gat2$parameter, levels=c("X2","X1"), labels = c("sigma[L]^2", "sigma[H]^2"))}
# if(K == 3){  gat2$parameter<-factor(gat2$parameter, levels=c("X3", "X2","X1"), labels = c("sigma[L]^2","sigma[M]^2","sigma[H]^2"))}
# 
# gat1$parameter<-factor(gat1$parameter)
# gat1$parameter<-factor(gat1$parameter, labels = c("gamma[0]", "gamma[1]", "phi", "tau"))
# 
# gat<-rbind(gat1, gat2)
# gat<-gat[!gat$parameter %in% c("alpha", "beta" , "tau"), ]
# 
# make_label <- function(value) {
#   x <- as.character(value)
#   bquote(italic(.(x))~subjects)
# }
# 
# plot_labeller <- function(variable, value) {
#   do.call(expression, lapply(levels(value), make_label))
# }
# 
# 
# gat<-gat[,1:2]
# 
# colnames(gat)[1]<-"level"
# 
# gat$type <- "posterior"
# if(K == 2){gat2$parameter<-factor(gat2$parameter, levels=c("X2","X1"), labels = c("sigma[L]^2",  "sigma[H]^2"))}
# if(K == 3){gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[M]^2", "sigma[H]^2"))}
# 
# gat_it<- gat
# gat_it$type <- "posterior"
# gat_it$country <- "Italy"
# 
# result7<- result$xi_ite
# result7<-lapply(result7, as.data.frame)
# result7<- data.frame(rbindlist(result7))
# result7$t <- 1:Time
# result7$ite<- rep(1:45000, each = Time)
# 
# # colnames(result7)<-c("state1", "state2", "state3", "t", "ite")
# colnames(result7)<-c(paste0("state",1:K),  "t", "ite")
# 
# result7<-data.frame(result7)
# 
# xi_it_sub<-result7[result7$ite %in% seq(35000, 45000, 10),]
# xi_it_agg<-aggregate(xi_it_sub[,1:K], by = list(t = xi_it_sub$t), FUN = mean)
# xi_it_agg$State<-  apply(xi_it_agg[,(2:(K+1))],1, which.max)
# # xi_it_agg$State<- factor(xi_it_agg$state1)
# # xi_it_agg$State<- factor(xi_it_agg$State, labels = c("L", "H"))
# xi_it_agg$t <- dates
# 
# if(K == 2){xi_it_agg$State<-factor((K+1)- xi_it_agg$State, labels = c("L", "H"))}
# if(K == 3){xi_it_agg$State<-factor((K+1)-xi_it_agg$State, labels = c("L","M", "H"))}
# 
# q<- ggplot(xi_it_agg)+labs(x ="time", y = "State", colour = "State")
# q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
# q<- q + geom_rect(data = xi_it_agg , aes(xmin = t, xmax = lead(t),  ymax  = State),  ymin= "L", fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
# # q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
# #q<- q + geom_step(aes(x = t, y = State, group = 1), color ="black", shape = 15 , size = 0.2)
# 
# # q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
# #                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# # 
# q<- q + labs(title = "Italy",
#              x = "Time", y = "State")+ theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                    size = 36, hjust = 0.5),
#                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                
#                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                axis.text=element_text(size=16),
#                                                                axis.line = element_blank(),
#                                                                axis.ticks = element_line(),
#                                                                panel.grid = element_blank(),
#                                                                legend.title = element_text(face="italic", size=16),
#                                                                legend.text = element_text(size=14),
#                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# q_it<- q
# 
# 
# result6<- data.frame(result$P_ite)
# result6<- result6[seq(35000, 45000, 10),1:(K*K)]
# colnames(result6)<-paste("p",rep(1:K, each = length(1:K)), 1:K, sep = "")
# ps<- apply(result6,2, FUN = mean)
# ps
# ps_m<-matrix(ps, K,K, byrow = T)
# ps_m<-round(ps_m,3)
# 
# 
# if(K == 2){
#   
#   rownames(ps_m)<- c("H", "L")
#   colnames(ps_m)<- c("H", "L")
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L", "H"))
#   
# }
# 
# if(K == 3){
#   
#   rownames(ps_m)<- c("H", "M", "L")
#   colnames(ps_m)<- c("H", "M", "L")
#   
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H","M", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L","M", "H"))
# }
# 
# # Create a ggheatmap
# ggheatmap1 <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
#   geom_tile(color = "white")+ labs(x= "States", y ="States", title = "Italy" )+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   theme_minimal() + 
#   geom_text(aes(Var1, Var2, label = value), color = "black", size = 8) + theme(legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                                  size = 22, hjust = 0.5),
#                                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                                
#                                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                                axis.text=element_text(size=16),
#                                                                                axis.line = element_blank(),
#                                                                                axis.ticks = element_line(),
#                                                                                panel.grid = element_blank(),
#                                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# ggheatmap_it <- ggheatmap1
# 
# EY = rep(0, Time)
# Var = rep(0, Time)
# DI = rep(0, Time)
# 
# agg<-aggregate(EL_x$w, by = list( i = EL_x$i, t = EL_x$t), FUN = sum) # strength distr through time
# 
# agg3<-aggregate(agg$x, by = list(agg$t), FUN = mean)
# agg4<-aggregate(agg$x, by = list(agg$t), FUN = var)
# 
# EY  = agg3$x
# Var = agg4$x
# DI = agg4$x/agg3$x
# 
# lkMetrics<- logPredictiveScoreMS_metrics(EL_princ$w,  result$beta_it, result$zi_ite,   result$xi_ite , Npages,  Time, K, D, 35*10^3, 45*10^3 , 10)
# 
# ResList_plot_db_long_a<-gather(data.frame(lkMetrics$EY), key = "Time", value = "Value")
# ResList_plot_db_long_a$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_a$Metric<- "Expected Strength"
# 
# ResList_plot_db_long_b<-gather(data.frame(lkMetrics$Var), key = "Time", value = "Value")
# ResList_plot_db_long_b$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_b$Metric<- "Variance Strength"
# 
# ResList_plot_db_long_c<-gather(data.frame(lkMetrics$DI), key = "Time", value = "Value")
# ResList_plot_db_long_c$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_c$Metric<- "Dispersion Index"
# 
# ResList_plot_db_long<-rbind(ResList_plot_db_long_a, ResList_plot_db_long_b, ResList_plot_db_long_c)
# ResList_plot_db_long$Metric<-factor(ResList_plot_db_long$Metric, levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ) )
# 
# q<-ggplot(ResList_plot_db_long, aes(x = Time, y = Value, group = Time), col ="red") + facet_wrap( vars(Metric), scales = "free")
# #p <- p +  geom_step(data = line_db, aes(x = date, y = value))
# q <- q +  geom_line(data = data.frame(Time = dates , Value = EY ,  Metric =   factor("Expected Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) ,col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = Var ,  Metric =   factor("Variance Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) , col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = DI ,  Metric =   factor("Dispersion Index", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value), col ="grey",  inherit.aes = F)
# q <- q +  geom_boxplot( fill = "red", col ="darkred", outlier.shape  = NA, lwd = 0.1)
# # p <- p +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
# #                         geom="crossbar", width=0.5)
# q<- q+ labs(y = "", title = "Italy") + theme_bw()
# q<- q + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                   size = 22, hjust = 0.5),
#               strip.text.x = element_text(size = 16,face = "italic"),
#               
#               axis.title = element_text( face = "italic",size = 16),
#               axis.title.y = element_text(size = 18, angle=90,vjust =2),
#               axis.title.x = element_text(size = 18, vjust = -0.2),
#               axis.text=element_text(size=8),
#               axis.line = element_blank(),
#               axis.ticks = element_line(),
#               panel.grid = element_blank(),
#               legend.title = element_text(face="italic", size=16),
#               legend.text = element_text(size=14),
#               panel.border = element_rect(colour = "black", fill = NA))
# gv_it<-q
# 
# interval = 25000:45000
# 
# ########CHANGE YOUR PATH ###########
# setwd("~/Repository/")
# #####################################
# 
# load("Data/Dynamic/DataEnv_SP_all.RData")
#
# ######## LOAD RESULTS  ############
# # load results.Rdata from Dynamic_01_Results_SP_extended.R
# load("")
# ###################################
# 
# sourceCpp("Code/Model/Extended_MS_LS_FE.cpp")
# sourceCpp("Code/Model/Extended_Predictive.cpp")
# 
# result$DIC
# a<-data.frame(name= unique(EL_x$i), max_czeros = 0  )
# 
# 
# for(i in 1:length(unique(EL_x$i) )){
#   
#   name = unique(EL_x$i)[i]
#   resa<-aggregate(EL_x$w[EL_x$i == name ], by = list(EL_x$t[EL_x$i == name ]), sum)
#   x<- rle(resa$x==0)
#   
#   if(length(x$lengths[x$values == TRUE])>0){
#     a[i,2]= max(x$lengths[x$values == TRUE])}
#   print(i)
# }
# 
# 
# thr <- c(a$name[a$max_czeros > 15])
# 
# 
# DBplane<-DBplane[! DBplane$fb_name %in% thr ,  ]
# DBplane$i<- as.numeric(factor(DBplane$fb_name))
# #
# EL_x<-EL_x[! EL_x$i %in% thr ,  ]
# EL_x<-EL_x[! EL_x$j %in% thr ,  ]
# #
# EL_princ<-EL_princ[! EL_princ$i %in% thr ,  ]
# EL_princ<-EL_princ[! EL_princ$j %in% thr ,  ]
# #
# EL<-EL[! EL$i %in% thr ,  ]
# EL<-EL[! EL$j %in% thr ,  ]
# 
# EL$ith<-as.numeric(factor(EL$i) , levels = unique(EL$i))
# EL$jth<-as.numeric(factor(EL$j) , levels = unique(EL$i))
# 
# EL_x$ith<-as.numeric(factor(EL_x$i) , levels = unique(EL_x$i))
# EL_x$jth<-as.numeric(factor(EL_x$j) , levels = unique(EL_x$i))
# 
# EL_princ$ith<-as.numeric(factor(EL_princ$i , levels = unique(EL_x$i)))
# EL_princ$jth<-as.numeric(factor(EL_princ$j, levels = unique(EL_x$i)))
# 
# pages_names_sp<- c(unique(EL_x$i)) # a vector of page names
# 
# #Main Dimensions
# N<- length(unique(EL_x$i_j)) #number of edges
# M<- length(unique(EL_x$i_j))/2 #number of unique edges
# Time<- length(unique(EL_x$t)) #number of unique dates
# Npages<- length(pages_names_sp) #number of unique pages
# K <- 2 #number of states
# D<-2
# 
# dates<- unique(as_date(EL_princ$tindex))
# 
# 
# top_names<- c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" )
# 
# beta_sp<- result$beta_it
# agg_beta_sp<-colMeans(beta_sp[interval, ])
# # 
# #################
# k = 1
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# #################
# 
# k = 1
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_1 = rbindlist(zi_res)
# zi_res_long_1 = data.frame(zi_res_long_1)
# zi_res_long_1$i = 1:Npages
# zi_res_long_1$it<- rep(1:length(zi_res), each = Npages)
# 
# agg_zeta_sp_a= aggregate(zi_res_long_1[,1:D], by = list(i= zi_res_long_1$i), FUN = mean)
# agg_zeta_sp_a$names = pages_names_sp
# colnames(agg_zeta_sp_a)[2:(D+1)]<-paste0("X",1:D)
# 
# #################
# k = 2
# subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
# subz_k_ref = as.matrix(subz_k[[length(interval)]])
# 
# zi_res<-subz_k
# zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
# zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
# 
# zi_res_long_2 = rbindlist(zi_res)
# zi_res_long_2 = data.frame(zi_res_long_2)
# zi_res_long_2$i = 1:Npages
# 
# agg_zeta_sp_b= aggregate(zi_res_long_2[,1:D], by = list(i= zi_res_long_2$i), FUN = mean)
# agg_zeta_sp_b$names = pages_names_sp
# colnames(agg_zeta_sp_b)[2:(D+1)]<-paste0("X",1:D)
# 
# 
# ############
# 
# if(K == 3){
#   k = 3
#   subz_k =  unlist(lapply(result$zi_ite[interval],"[", k), recursive = F)
#   subz_k_ref = as.matrix(subz_k[[length(interval)]])
#   
#   zi_res<-subz_k
#   zi_res<-lapply(subz_k, FUN = function(x){procrustes_cpp(as.matrix(subz_k_ref),as.matrix(x) )}) #apply procrustes trasnformation
#   zi_res<-lapply(zi_res, FUN = function(x){as.data.frame(x)}) #apply procrustes trasnformation
#   
#   zi_res_long_3 = rbindlist(zi_res)
#   zi_res_long_3 = data.frame(zi_res_long_3)
#   zi_res_long_3$i = 1:Npages
#   
#   agg_zeta_sp_c= aggregate(zi_res_long_3[,1:D], by = list(i= zi_res_long_3$i), FUN = mean)
#   agg_zeta_sp_c$names = pages_names_sp
#   colnames(agg_zeta_sp_c)[2:(D+1)]<-paste0("X",1:D)
#   
#   
# }
# 
# ##############
# 
# pew_sp_a<- agg_zeta_sp_a[agg_zeta_sp_a$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" ), ]
# pew_sp_a$State<- "State H"
# pew_sp_a$score_a<- c(4.5, 3.8, 4.2, 3.4, 3.5) - mean( c(4.5, 3.8, 4.2, 3.4, 3.5))
# pew_sp_a$score_b<- c(3.3, 3.1, 3.2, 3.1, 2.8) - mean( c(3.3, 3.1, 3.2, 3.1, 2.8) )
# 
# cor(pew_sp_a$X1, pew_sp_a$score_a)
# 
# pew_sp_b<- agg_zeta_sp_b[agg_zeta_sp_b$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" ), ]
# pew_sp_b$State<- "State L"
# pew_sp_b$score_a<- c(4.5, 3.8, 4.2, 3.4, 3.5) - mean( c(4.5, 3.8, 4.2, 3.4, 3.5))
# pew_sp_b$score_b<- c(3.3, 3.1, 3.2, 3.1, 2.8) - mean( c(3.3, 3.1, 3.2, 3.1, 2.8) )
# pew_sp<-rbind(pew_sp_a, pew_sp_b)
# 
# if(K == 3){
#   pew_sp_b$State<- "State M"
#   pew_sp_c<- agg_zeta_sp_c[agg_zeta_sp_c$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" ), ]
#   pew_sp_c$State<- "State L"
#   pew_sp_c$score_a<- c(4.5, 3.8, 4.2, 3.4, 3.5) - mean( c(4.5, 3.8, 4.2, 3.4, 3.5))
#   pew_sp_c$score_b<- c(3.3, 3.1, 3.2, 3.1, 2.8) - mean( c(3.3, 3.1, 3.2, 3.1, 2.8) )
#   pew_sp<-rbind(pew_sp_a, pew_sp_b, pew_sp_c)
# }
# 
# pew_sp$country<- "Spain"
# 
# db_sp_a<-data.frame(agg_zeta_sp_a[,2:(D+1)])
# colnames(db_sp_a)   <- paste0("agg_zeta_sp_", 1:D)             
# db_sp_a$agg_beta_sp= agg_beta_sp
# db_sp_a$names=pages_names_sp
# 
# db_sp_a$State<- "State H"
# 
# db_sp_b<-data.frame(agg_zeta_sp_b[,2:(D+1)])
# colnames(db_sp_b)   <- paste0("agg_zeta_sp_", 1:D)             
# db_sp_b$agg_beta_sp= agg_beta_sp
# db_sp_b$names=pages_names_sp
# 
# db_sp_b$State<- "State L"
# 
# db_sp<- rbind(db_sp_a,db_sp_b)
# 
# 
# if(K ==3){
#   db_sp_b$State<- "State M"
#   
#   db_sp_c<-data.frame(agg_zeta_sp_c[,2:(D+1)])
#   colnames(db_sp_c)   <- paste0("agg_zeta_sp_", 1:D)             
#   db_sp_c$agg_beta_sp= agg_beta_sp
#   db_sp_c$names=pages_names_sp
#   
#   db_sp_c$State<- "State L"
#   db_sp<- rbind(db_sp_a,db_sp_b, db_sp_c)
# }
# 
# if(K == 2){db_sp$State<-factor(db_sp$State,  levels = c("State L", "State H"))}
# if(K == 3){db_sp$State<-factor(db_sp$State,  levels = c("State L", "State M", "State H"))}
# 
# p2<-ggplot(db_sp, aes(x = agg_zeta_sp_1, y = agg_zeta_sp_2,  label = names)) +facet_wrap(~State, scales = "free_y")
# p2<- p2+ geom_point(col = "#ff420f", aes( size = exp(agg_beta_sp),))
# p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
# #p2<- p2 + geom_text_repel(data = db_sp, size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Latent Coordinate", title = "Spain", col = "white")
# #p2<- p2 + geom_text_repel(data = db_sp[db_sp$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "Spain", col = "white", size = "Individual Effect")
# p2<- p2 + geom_text_repel(data = db_sp,size = 1,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "2nd Latent Coordinate", title = "Spain", col = "white", size = "Individual Effect")
# p2 <- p2 + theme(panel.grid = element_blank())  
# p2 <- p2 + theme_minimal() + theme(strip.placement = "outside",legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                  size = 14),
#                                    strip.text.x = element_blank(),
#                                    
#                                    axis.title = element_text( face = "italic",size = rel(1)),
#                                    axis.title.y = element_text(size = 14, angle=90,vjust =2),
#                                    axis.title.x = element_text(size = 14, vjust = -0.2),
#                                    axis.text=element_text(size=12),
#                                    axis.line = element_blank(),
#                                    axis.ticks = element_line(),
#                                    panel.grid = element_blank(),
#                                    panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# 
# 
# p2_sp<-p2
# 
# result3<- result$phi_gamma_ite
# 
# colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
# result3<-data.frame(result3)
# 
# 
# gat1<-gather(result3[result3$ite %in% interval,1:4], key = "parameter", value = "value")
# 
# result8<-data.frame(result$sigma_ite)
# result8$X1<- result8$X1^2
# result8$X2<- result8$X2^2
# if(K == 3){result8$X3<- result8$X3^2}
# 
# result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
# result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
# if(K == 3){result8$sc_mean <- cumsum(result8$X3)/(1:length(result8$X3))}
# 
# result8$ite<- 1:length(result8$X2)
# 
# 
# gat2<-gather(result8[result8$ite %in% interval,1:K], key = "parameter", value = "value")
# if(K == 2){gat2$parameter<-factor(gat2$parameter, levels=c("X2","X1"), labels = c("sigma[L]^2", "sigma[H]^2"))}
# if(K == 3){gat2$parameter<-factor(gat2$parameter, levels=c("X3", "X2","X1"), labels = c("sigma[L]^2","sigma[M]^2","sigma[H]^2"))}
# 
# gat1$parameter<-factor(gat1$parameter)
# gat1$parameter<-factor(gat1$parameter, labels = c("gamma[0]", "gamma[1]", "phi", "tau"))
# 
# gat<-rbind(gat1, gat2)
# gat<-gat[!gat$parameter %in% c("alpha", "beta" , "tau"), ]
# 
# make_label <- function(value) {
#   x <- as.character(value)
#   bquote(italic(.(x))~subjects)
# }
# 
# plot_labeller <- function(variable, value) {
#   do.call(expression, lapply(levels(value), make_label))
# }
# 
# 
# gat<-gat[,1:2]
# # gat_prior<-gat
# # gat_prior$value <-c(prior_phi, prior_g0, prior_g1)
# # gat_prior$type <- "prior"
# 
# colnames(gat)[1]<-"level"
# 
# gat$type <- "posterior"
# 
# if(K == 2){gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[H]^2"))}
# if(K == 3){gat$level<-factor(gat$level, levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[M]^2", "sigma[H]^2"))}
# 
# gat_sp<- gat
# gat_sp$country <- "Spain"
# 
# 
# 
# result7<- result$xi_ite
# result7<-lapply(result7, as.data.frame)
# result7<- data.frame(rbindlist(result7))
# result7$t <- 1:Time
# result7$ite<- rep(1:45000, each = Time)
# 
# colnames(result7)<-c(paste0("state",1:K), "t", "ite")
# result7<-data.frame(result7)
# 
# xi_sp_sub<-result7[result7$ite %in% seq(35000, 45000, 10),]
# xi_sp_agg<-aggregate(xi_sp_sub[,1:K], by = list(t = xi_sp_sub$t), FUN = mean)
# xi_sp_agg$State<-  apply(xi_sp_agg[,2:(K+1)],1, which.max)
# 
# # xi_it_agg$State<- factor(xi_it_agg$state1)
# # xi_it_agg$State<- factor(xi_it_agg$State, labels = c("L", "H"))
# xi_sp_agg$t <- dates
# 
# if(K == 2){ xi_sp_agg$State<-factor(4-xi_sp_agg$State, labels = c("L", "H"))}
# if(K == 3){ xi_sp_agg$State<-factor(4-xi_sp_agg$State, labels = c("L","M", "H"))}
# 
# q<- ggplot(xi_sp_agg)+labs(x ="time", y = "State", colour = "State")
# q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
# # q<- q + geom_rect(data = xi_it_agg_area , aes(xmin = t, xmax = lead,  ymax  = state1),  ymin= min(as.numeric(xi_it_agg$state1)), fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
# # q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
# #q<- q + geom_step(aes(x = t, y = (State), group = 1), color ="black", shape = 15 , size = 0.2)
# q<- q + geom_rect(data = xi_sp_agg , aes(xmin = t, xmax = lead(t),  ymax  = State),  ymin= "L", fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
# # q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
# #                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# # 
# q<- q + labs(title = "Spain",
#              x = "Time", y = "State")+ theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                    size = 36, hjust = 0.5),
#                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                
#                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                axis.text=element_text(size=16),
#                                                                axis.line = element_blank(),
#                                                                axis.ticks = element_line(),
#                                                                panel.grid = element_blank(),
#                                                                legend.title = element_text(face="italic", size=16),
#                                                                legend.text = element_text(size=14),
#                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# q_sp<- q
# 
# result6<- data.frame(result$P_ite)
# result6<- result6[seq(35000, 45000, 10),1:(K*K)]
# colnames(result6)<-paste("p",rep(1:K, each = length(1:K)), 1:K, sep = "")
# ps<- apply(result6,2, FUN = mean)
# ps
# ps_m<-matrix(ps, K,K, byrow = T)
# ps_m<-round(ps_m,3)
# 
# 
# if(K == 2){
#   
#   rownames(ps_m)<- c("H", "L")
#   colnames(ps_m)<- c("H", "L")
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L", "H"))
#   
# }
# 
# if(K == 3){
#   
#   rownames(ps_m)<- c("H", "M", "L")
#   colnames(ps_m)<- c("H", "M", "L")
#   
#   
#   melted_cormat <- reshape2::melt(t(ps_m), na.rm = TRUE)
#   melted_cormat$value<-round(melted_cormat$value,3)
#   melted_cormat$Var1<-factor(melted_cormat$Var1, levels = c("H","M", "L"))
#   melted_cormat$Var2<-factor(melted_cormat$Var2, levels = c("L","M", "H"))
# }
# 
# # Create a ggheatmap
# ggheatmap1 <- ggplot(melted_cormat, aes(Var1, Var2, fill = value))+
#   geom_tile(color = "white")+ labs(x= "States", y ="States", title = "Spain" )+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0, limit = c(-1,1), space = "Lab", 
#                        name="Pearson\nCorrelation") +
#   theme_minimal() + 
#   geom_text(aes(Var1, Var2, label = value), color = "black", size = 8) + theme(legend.position="none", plot.title = element_text(face = "italic",
#                                                                                                                                  size = 22, hjust = 0.5),
#                                                                                strip.text.x = element_text(size = 24,face = "italic"),
#                                                                                
#                                                                                axis.title = element_text( face = "italic",size = rel(1)),
#                                                                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
#                                                                                axis.title.x = element_text(size = 22, vjust = -0.2),
#                                                                                axis.text=element_text(size=16),
#                                                                                axis.line = element_blank(),
#                                                                                axis.ticks = element_line(),
#                                                                                panel.grid = element_blank(),
#                                                                                panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# ggheatmap_sp <- ggheatmap1
# 
# EY = rep(0, Time)
# Var = rep(0, Time)
# DI = rep(0, Time)
# 
# agg<-aggregate(EL_x$w, by = list( i = EL_x$i, t = EL_x$t), FUN = sum) # strength distr through time
# 
# agg3<-aggregate(agg$x, by = list(agg$t), FUN = mean)
# agg4<-aggregate(agg$x, by = list(agg$t), FUN = var)
# 
# EY  = agg3$x
# Var = agg4$x
# DI = agg4$x/agg3$x
# 
# lkMetrics<- logPredictiveScoreMS_metrics(EL_princ$w,  result$beta_it, result$zi_ite,   result$xi_ite , Npages,  Time, K, D, 35*10^3, 45*10^3 , 10)
# 
# ResList_plot_db_long_a<-gather(data.frame(lkMetrics$EY), key = "Time", value = "Value")
# ResList_plot_db_long_a$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_a$Metric<- "Expected Strength"
# 
# ResList_plot_db_long_b<-gather(data.frame(lkMetrics$Var), key = "Time", value = "Value")
# ResList_plot_db_long_b$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_b$Metric<- "Variance Strength"
# 
# ResList_plot_db_long_c<-gather(data.frame(lkMetrics$DI), key = "Time", value = "Value")
# ResList_plot_db_long_c$Time<- rep(dates, each = length(seq(35*10^3, 45*10^3 , 10)))
# ResList_plot_db_long_c$Metric<- "Dispersion Index"
# 
# ResList_plot_db_long<-rbind(ResList_plot_db_long_a, ResList_plot_db_long_b, ResList_plot_db_long_c)
# ResList_plot_db_long$Metric<-factor(ResList_plot_db_long$Metric, levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ) )
# 
# q<-ggplot(ResList_plot_db_long, aes(x = Time, y = Value, group = Time), col ="red") + facet_wrap( vars(Metric), scales = "free")
# #p <- p +  geom_step(data = line_db, aes(x = date, y = value))
# q <- q +  geom_line(data = data.frame(Time = dates , Value = EY ,  Metric =   factor("Expected Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) ,col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = Var ,  Metric =   factor("Variance Strength", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value) , col ="grey", inherit.aes = F)
# q <- q +  geom_line(data = data.frame(Time = dates , Value = DI ,  Metric =   factor("Dispersion Index", levels = c( "Expected Strength", "Variance Strength", "Dispersion Index" ))), aes(x = Time, y = Value), col ="grey",  inherit.aes = F)
# q <- q +  geom_boxplot( fill = "red", col ="darkred", outlier.shape  = NA, lwd = 0.1)
# # p <- p +  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
# #                         geom="crossbar", width=0.5)
# q<- q+ labs(y = "", title = "Spain") + theme_bw()
# q<- q + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                   size = 22, hjust = 0.5),
#               strip.text.x = element_text(size = 16,face = "italic"),
#               
#               axis.title = element_text( face = "italic",size = 16),
#               axis.title.y = element_text(size = 18, angle=90,vjust =2),
#               axis.title.x = element_text(size = 18, vjust = -0.2),
#               axis.text=element_text(size=8),
#               axis.line = element_blank(),
#               axis.ticks = element_line(),
#               panel.grid = element_blank(),
#               legend.title = element_text(face="italic", size=16),
#               legend.text = element_text(size=14),
#               panel.border = element_rect(colour = "black", fill = NA))
# gv_sp<-q
# 
# 
# ####
# PEW_all<-rbind(pew_de, pew_fr, pew_it, pew_sp)
# if(K == 3){PEW_all$State<- factor(PEW_all$State, levels = c("State L", "State M", "State H")  )}
# if(K == 2){PEW_all$State<- factor(PEW_all$State, levels = c("State L", "State H")  )}
# 
# ##########
# #Figure 8 
# 
# gat_all<-rbind(gat_fr, gat_de, gat_it, gat_sp)
# 
# sturges <- function(x){ pretty(range(x),
#                                n = nclass.Sturges(x),
#                                min.n = 1)}
# 
# gat_all<- gat_all[!gat_all$level %in% "gamma[1]" ,]
# gat_all$it<-25000:45000
# 
# z <- ggplot(gat_all[gat_all$it %in% seq(25000,45000,5),],  aes(x=value, group = country, linetype = country )) +  geom_density(fill = "#ff420f", alpha = 0.3,  adjust = 3.5)  + facet_wrap(.~level, scales = "free",  labeller = label_parsed)
# z<- z + geom_hline(yintercept = 0 , col = "white", linewidth = 1)
# # z <-  z+ geom_histogram(data = gat[gat$level == "gamma[0]",], aes( y = after_stat(density)), binwidth=0.1, fill = "#ff420f", alpha = 0.6 , col ="black", breaks = sturges(gat[gat$level == "gamma[0]","value"] ) ) 
# # z <- z+ geom_histogram(data = gat[gat$level == "gamma[1]",],  aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "gamma[1]","value"] ) ) 
# # z <- z+ geom_histogram(data = gat[gat$level == "phi",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "phi","value"] ) ) 
# # z <- z+ geom_histogram(data = gat[gat$level == "sigma[L]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[L]^2","value"] ) ) 
# # z <- z+ geom_histogram(data = gat[gat$level == "sigma[H]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[H]^2","value"] ) ) 
# # z<- z + geom_density(linetype = 1,  n = 10000, adjust = 2)
# z <- z + scale_linetype_manual(values=c("solid", "twodash", "dotted", "dashed")) 
# z<-  z + stat_function( data = data.frame(value = 0, yaxis = 0, type = "prior", level = factor("phi",levels = levels(gat_all$level)) ), fun = dgamma, args = list(shape=0.01, rate =0.01), geom = "area", fill ="black", alpha = 0.3, linetype = "dashed", colour = "red", inherit.aes = F)
# z<-  z + stat_function( data = data.frame(value = 0, yaxis = 0, type = "prior", level = factor("gamma[0]",levels = levels(gat_all$level)
# ) ),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3, linetype = "dashed", colour = "red", inherit.aes = F)
# # z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, type = "prior", level = factor("gamma[1]",levels = levels(gat_all$level))),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, type = "prior", level = factor("sigma[L]^2",levels = levels(gat_all$level))),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)
# if(K == 3){z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, type = "prior", level = factor("sigma[M]^2",levels = levels(gat_all$level))),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)}
# z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, type = "prior", level = factor("sigma[H]^2",levels = levels(gat_all$level))),fun = dinvgamma, args = list(shape = 0.01, scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)
# #z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed)
# # z<- z + scale_shape_manual(labels =  parse_format())
# z <- z + labs(x = "Value",  y = "", linetype="Country")+  theme_minimal() + theme( legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "italic",
#                                                                                                                                                                         size = 12, hjust = 0.5),
#                                                                                    legend.key.size = unit(1, 'cm'),
#                                                                                    legend.title = element_text(size=20), #change legend title font size
#                                                                                    legend.text = element_text(size=18),
#                                                                                    axis.text=element_text(size=8),
#                                                                                    strip.text.x = element_text(size = 18, face = "italic"),
#                                                                                    axis.title.x = element_blank(),
#                                                                                    axis.title.y = element_text(size = 22, angle=90,vjust =2, face = "italic"),
#                                                                                    axis.line = element_blank(),
#                                                                                    axis.ticks = element_line(),
#                                                                                    axis.ticks.x  = element_line(),
#                                                                                    panel.spacing.x = unit(4, "mm"),
#                                                                                    panel.grid = element_blank(),
#                                                                                    panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# z<-z + guides(linetype=guide_legend(ncol=4))
# 
# w_all<-z
# 
# ggsave(w_all, filename = paste0("Figures/Dynamic/Figure8_v2_extended_gamma_K_", K,".pdf"), units = "cm", width = 16*2, height = 9*2 )
# 
# shift_legend <- function(p){
#   
#   # check if p is a valid object
#   if(!"gtable" %in% class(p)){
#     if("ggplot" %in% class(p)){
#       gp <- ggplotGrob(p) # convert to grob
#     } else {
#       message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
#       return(p)
#     }
#   } else {
#     gp <- p
#   }
#   
#   # check for unfilled facet panels
#   facet.panels <- grep("^panel", gp[["layout"]][["name"]])
#   empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
#   empty.facet.panels <- facet.panels[empty.facet.panels]
#   if(length(empty.facet.panels) == 0){
#     message("There are no unfilled facet panels to shift legend into. Returning original plot.")
#     return(p)
#   }
#   
#   # establish extent of unfilled facet panels (including any axis cells in between)
#   empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
#   empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
#                              max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
#   names(empty.facet.panels) <- c("t", "l", "b", "r")
#   
#   # extract legend & copy over to location of unfilled facet panels
#   guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
#   if(length(guide.grob) == 0){
#     message("There is no legend present. Returning original plot.")
#     return(p)
#   }
#   gp <- gtable_add_grob(x = gp,
#                         grobs = gp[["grobs"]][[guide.grob]],
#                         t = empty.facet.panels[["t"]],
#                         l = empty.facet.panels[["l"]],
#                         b = empty.facet.panels[["b"]],
#                         r = empty.facet.panels[["r"]],
#                         name = "new-guide-box")
#   
#   # squash the original guide box's row / column (whichever applicable)
#   # & empty its cell
#   guide.grob <- gp[["layout"]][guide.grob, ]
#   if(guide.grob[["l"]] == guide.grob[["r"]]){
#     gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
#   }
#   if(guide.grob[["t"]] == guide.grob[["b"]]){
#     gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
#   }
#   gp <- gtable_remove_grobs(gp, "guide-box")
#   
#   return(gp)
# }
# 
# 
# 
# 
# ## Initiate writing to PDF file
# pdf(paste0("Figures/Dynamic/Figure8_v2_extended_K", K,".pdf"),  width = 16, height = 9)
# 
# ## Create a graphical object g here
# grid.draw(shift_legend(w_all))
# # print it
# 
# ## Stop writing to the PDF file
# dev.off()
# 
# 
# ##########
# 
# p2_fr<-p2_fr+labs( x = "1st Coordinate (Leaning)", y= "2nd Coordinate") 
# 
# gp <- ggplotGrob(p2_fr)
# 
# # gp$layout #helps you to understand the gtable object 
# # gtable_show_layout(ggplotGrob(p)) #helps you to understand the gtable object 
# 
# t_title <- gp$layout[gp$layout[['name']] == 'title' ,][['t']]
# t_strip <- gp$layout[grepl('strip', gp$layout[['name']]),][['t']]
# 
# gp$layout[gp$layout[['name']] == 'title' ,][['t']] <- unique(t_strip)
# gp$layout[gp$layout[['name']] == 'title' ,][['b']] <- unique(t_strip)
# gp$layout[grepl('strip', gp$layout[['name']]),][['t']] <- t_title
# gp$layout[grepl('strip', gp$layout[['name']]),][['b']] <- t_title
# 
# gp_fr<-gp
# 
# 
# 
# p2_de<-p2_de+labs( x = "1st Coordinate (Leaning)", y= "2nd Coordinate") 
# 
# gp <- ggplotGrob(p2_de)
# 
# gp_de<-gp
# 
# 
# p2_it<-p2_it+labs( x = "1st Coordinate (Leaning)", y= "2nd Coordinate") 
# 
# 
# gp <- ggplotGrob(p2_it)
# 
# gp_it<-gp
# 
# 
# p2_sp<-p2_sp+labs( x = "1st Coordinate (Leaning)", y= "2nd Coordinate") 
# 
# 
# gp <- ggplotGrob(p2_sp)
# 
# gp_sp<-gp
# 
# 
# pp<- (p2_fr/p2_de/p2_it/p2_sp) & theme(axis.title = element_text(size = 16) ) & labs( x = "1st Coordinate (Leaning)", "2nd Coordinate")
# pp<- ((wrap_elements(full =  gp_fr))/(wrap_elements(full =  gp_de))/(wrap_elements(full =  gp_it))/(wrap_elements(full =  gp_sp)))  & theme(axis.title = element_text(size = 16) )
# 
# 
# ggsave(pp, filename = paste0("Figures/Dynamic/Figure9_extended_K",K,"_try.pdf"), units = "cm", width = 21, height = 29.7 )
# 
# 
# p5<- ggplot(PEW_all) + facet_wrap(.~State) +geom_point(aes(x =  X1 , y = score_a ), col = "#ff420f") +geom_text_repel(aes(x =  X1 , y = score_a, label = names ), size = 3) + theme_minimal()
# p5<- p5 + labs(x = "Latent Leaning", y = "PEW score") # + geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
# p5<- p5 + stat_cor(aes(x =  X1 , y = score_a, label = after_stat(r.label))) #+ xlim(-0.25, 0.25)  
# p5 <- p5 + theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                        size = 16),
#                                    strip.text.x = element_text(size = 16,face = "italic"),
#                                    
#                                    axis.title = element_text( face = "italic",size = rel(1)),
#                                    axis.title.y = element_text(size = 16, angle=90,vjust =2),
#                                    axis.title.x = element_text(size = 16, vjust = -0.2),
#                                    axis.text=element_text(size=14),
#                                    axis.line = element_blank(),
#                                    axis.ticks = element_line(),
#                                    panel.grid = element_blank(),
#                                    legend.title = element_text(face="italic", size=16),
#                                    legend.text = element_text(size=14),
#                                    panel.border = element_rect(colour = "black", fill = NA))
# 
# 
# 
# 
# 
# p5
# 
# 
# ggsave(p5, filename = paste0("Figures/Dynamic/Figure10_extended_K",K,".pdf"), units = "cm", width = 16*2, height = 9*2 )
# 
# 
# q_all <- (q_fr|q_de)/(q_it|q_sp)+(ggheatmap_fr|ggheatmap_de|ggheatmap_it|ggheatmap_sp) + plot_annotation(tag_levels = 'A') & theme(plot.title = element_text(size = 24), axis.text = element_text(size = 18)  )
# ggsave(q_all, filename = paste0("Figures/Dynamic/Figure11_extended_K",K,".pdf"), units = "cm", width = 16*3, height = 9*3 )
# 
# load("gv_fr.Rdata")
# load("gv_it.Rdata")
# load("gv_de.Rdata")
# load("gv_sp.Rdata")
# 
# gv<- (gv_fr|gv_de)/(gv_it|gv_sp)
# ggsave(gv, filename = paste0("Figures/Dynamic/Figure12_extended_K",K,".pdf"), units = "cm", width = 16*3, height = 9*3 )
# 
