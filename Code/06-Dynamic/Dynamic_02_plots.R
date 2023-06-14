
rm(list = ls())

#load the libraries

list.of.packages <- c("ggplot2", "patchwork", "ggrepel", "tidyr", "dplyr", "data.table", "patchwork", "scales", "lubridate", "ggpmisc", "MCMCpack", "tidyverse", "gtable","grid", "reshape2", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(patchwork)
library(ggrepel)
library(tidyr)
library(dplyr)
library(data.table)
library(patchwork)
library(scales)
library(lubridate)
library(ggpmisc)
library(MCMCpack)
library(tidyverse)
library(gtable)
library(grid)
library(reshape2)
library(ggpubr)



########CHANGE YOUR PATH ###########
setwd("~/Documents/GitHub/Dyn-MS-LS-Media/")
#####################################

source("Code/05-Static/Misc/CoordCartesian.R")


################GERMANY########################

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
#p2<- p2 + geom_text_repel(data = db_de,size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Germany", col = "white")
p2<- p2 + geom_text_repel(data = db_de[db_de$names %in% top_names,] ,size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Germany", col = "white")
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


gat1<-gather(result3[result3$ite %in% seq(25000, 35000, 10 ),1:4], key = "parameter", value = "value")

result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)


gat2<-gather(result8[result8$ite %in% seq(25000, 35000, 1 ),1:2], key = "parameter", value = "value")
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

gat_de<- gat
gat_de$country <- "Germany"


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



###################FRANCE###################

load("Data/Dynamic/DataEnv_FR_all.RData")
load("Code/06-Dynamic/Results/RESULT_SC_FRz_fe.RData")

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

pages_names_fr<- c(unique(EL_x$i)) # a vector of page names

dates<- unique(as_date(DBplane$t))

top_names<-c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart")

beta_fr<- result[[1]]

zeta_fr_a<- result[[3]]
zeta_fr_a<-data.frame(zeta_fr_a)

zeta_fr_b<- result[[4]]
zeta_fr_b<-data.frame(zeta_fr_b)


colnames(zeta_fr_a) <-c("X1","X2", "i", "it")
colnames(zeta_fr_b) <-c("X1","X2", "i", "it")


zeta_fr_a_sub <-zeta_fr_a[zeta_fr_a$it %in% seq(20000, 35000, 1), ]
zeta_fr_b_sub <-zeta_fr_b[zeta_fr_b$it %in% seq(20000, 35000, 1), ]


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


agg_zeta_fr_a<-aggregate(zeta_fr_a_sub[,1:2], by = list(zeta_fr_a_sub$i), FUN = mean )
agg_zeta_fr_a$X1<-colMeans(zi_a_it[seq(20000, 35000, ), ])

agg_beta_fr<-colMeans(beta_fr[seq(20000, 35000, 1), ])
agg_zeta_fr_a$names<- pages_names_fr

agg_zeta_fr_b<-aggregate(zeta_fr_b_sub[,1:2], by = list(zeta_fr_b_sub$i), FUN = mean )
agg_zeta_fr_b$X1<-colMeans(zi_b_it[seq(20000, 35000, ), ])

agg_beta_fr<-colMeans(beta_fr[seq(20000, 35000, 1), ])
agg_zeta_fr_b$names<- pages_names_fr


pew_fr_a<- agg_zeta_fr_a[agg_zeta_fr_a$names %in%c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
pew_fr_a$State<- "State H"
pew_fr_a$score_a<-  c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1) - mean( c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1))
pew_fr_a$score_b<-  c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3) - mean( c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3) )

cor(pew_fr_a$X1, pew_fr_a$score_a)

pew_fr_b<- agg_zeta_fr_b[agg_zeta_fr_b$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
pew_fr_b$State<- "State L"
pew_fr_b$score_a<-  c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1) - mean(c(3.6, 3.2, 3.2, 2.7, 3))
pew_fr_b$score_b<- c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3) - mean(c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3))

pew_fr<-rbind(pew_fr_a, pew_fr_b)
cor(pew_fr_b$X1, pew_fr_b$score_a)
pew_fr$country<- "France"

#####################

db_fr_a<-data.frame(agg_zeta_fr= agg_zeta_fr_a$X1,agg_beta_fr= agg_beta_fr, names=pages_names_fr )
db_fr_a$State<- "State H"
db_fr_b<-data.frame(agg_zeta_fr= agg_zeta_fr_b$X1,agg_beta_fr= agg_beta_fr, names=pages_names_fr )
db_fr_b$State<- "State L"



db_fr<- rbind(db_fr_a,db_fr_b)
db_fr$State<-factor(db_fr$State,  levels = c("State L", "State H"))

#credible ellipse state a

p2<-ggplot(db_fr, aes(x = agg_zeta_fr, y = agg_beta_fr, label = names)) +facet_wrap(~State, scales = "free_y")
p2<- p2+ geom_point(col = "#ff420f", size = 1)  
p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
p2<- p2 + geom_text_repel(data = db_fr[db_fr$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "France", col = "white")
p2 <- p2 + theme(panel.grid = element_blank())  
p2 <- p2 + theme_minimal() + theme(strip.placement = "outside",legend.position="none", plot.title = element_text(face = "italic", size = 15),
                                   
                                   
                                   axis.title = element_text( face = "italic",size = rel(1)),
                                   axis.title.y = element_text(size = 14, angle=90,vjust =2),
                                   axis.title.x = element_blank(),
                                   axis.text=element_text(size=12),
                                   axis.line = element_blank(),
                                   axis.ticks = element_line(),
                                   panel.grid = element_blank(),
                                   panel.border = element_rect(colour = "black", fill = NA),
                                   strip.text.x = element_text(size = 14,face = "italic")  )



p2_fr<-p2


gp <- ggplotGrob(p2_fr)

# gp$layout #helps you to understand the gtable object 
# gtable_show_layout(ggplotGrob(p)) #helps you to understand the gtable object 

t_title <- gp$layout[gp$layout[['name']] == 'title' ,][['t']]
t_strip <- gp$layout[grepl('strip', gp$layout[['name']]),][['t']]

gp$layout[gp$layout[['name']] == 'title' ,][['t']] <- unique(t_strip)
gp$layout[gp$layout[['name']] == 'title' ,][['b']] <- unique(t_strip)
gp$layout[grepl('strip', gp$layout[['name']]),][['t']] <- t_title
gp$layout[grepl('strip', gp$layout[['name']]),][['b']] <- t_title



result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)


gat1<-gather(result3[result3$ite %in% seq(25000, 35000, 10 ),1:4], key = "parameter", value = "value")

result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)


gat2<-gather(result8[result8$ite %in% seq(25000, 35000, 1 ),1:2], key = "parameter", value = "value")
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

gat_fr<- gat
gat_fr$country <- "France"

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
z <- z + labs(title = "France", 
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




w_fr<-z 



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
xi_it_agg_area$lead[length(xi_it_agg_area$lead)]<- xi_it_agg$t[length(xi_it_agg$t)] +1

q<- ggplot(xi_it_agg)+labs(x ="time", y = "State", colour = "State")
#q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
q<- q + geom_rect(data = xi_it_agg_area , aes(xmin = t, xmax = lead,  ymax  = state1),  ymin= min(as.numeric(xi_it_agg$state1)), fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
#q<- q + geom_step(aes(x = t, y = as.numeric(State)), color ="black", shape = 15 , size = 0.2)

# q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
#                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# 
q<- q + labs(title = "France",
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



q_fr<- q





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
  geom_tile(color = "white")+ labs(x= "States", y ="States", title = "France" )+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ coord_flip()+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 8) + theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                                                                 size = 22, hjust = 0.5),
                                                                               strip.text.x = element_text(size = 24,face = "italic"),
                                                                               strip.placement = 'outside',
                                                                               
                                                                               axis.title = element_text( face = "italic",size = rel(1)),
                                                                               axis.title.y = element_text(size = 22, angle=90,vjust =2),
                                                                               axis.title.x = element_text(size = 22, vjust = -0.2),
                                                                               axis.text=element_text(size=16),
                                                                               axis.line = element_blank(),
                                                                               axis.ticks = element_line(),
                                                                               panel.grid = element_blank(),
                                                                               panel.border = element_rect(colour = "black", fill = NA))



ggheatmap_fr <- ggheatmap1


###########################


load("Data/Dynamic/DataEnv_IT_all.RData")
load("Code/06-Dynamic/Results/RESULT_SC_ITz_fe.RData")


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

pages_names_it<- c(unique(EL_x$i)) # a vector of page names


dates<- unique(as_date(DBplane$t))


top_names<-c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" )

beta_it<- result[[1]]

zeta_it_a<- result[[3]]
zeta_it_a<-data.frame(zeta_it_a)

zeta_it_b<- result[[4]]
zeta_it_b<-data.frame(zeta_it_b)


colnames(zeta_it_a) <-c("X1","X2", "i", "it")
colnames(zeta_it_b) <-c("X1","X2", "i", "it")


zeta_it_a_sub <-zeta_it_a[zeta_it_a$it %in% seq(20000, 35000, 1), ]
zeta_it_b_sub <-zeta_it_b[zeta_it_b$it %in% seq(20000, 35000, 1), ]


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

agg_zeta_it_a<-aggregate(zeta_it_a_sub[,1:2], by = list(zeta_it_a_sub$i), FUN = mean )
agg_zeta_it_a$X1<-colMeans(zi_a_it[seq(20000, 35000, ), ])

agg_beta_it<-colMeans(beta_it[seq(20000, 35000, 1), ])
agg_zeta_it_a$names<- pages_names_it

agg_zeta_it_b<-aggregate(zeta_it_b_sub[,1:2], by = list(zeta_it_b_sub$i), FUN = mean )
agg_zeta_it_b$X1<-colMeans(zi_b_it[seq(20000, 35000, ), ])

agg_beta_it<-colMeans(beta_it[seq(20000, 35000, 1), ])
agg_zeta_it_b$names<- pages_names_it


pew_it_a<- agg_zeta_it_a[agg_zeta_it_a$names %in%c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ), ]
pew_it_a$State<- "State H"
pew_it_a$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5) - mean( c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5))
pew_it_a$score_b<-  c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7) - mean( c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)  )

cor(pew_it_a$X1, pew_it_a$score_a)

pew_it_b<- agg_zeta_it_b[agg_zeta_it_b$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ), ]
pew_it_b$State<- "State L"
pew_it_b$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5) - mean( c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5))
pew_it_b$score_b<-  c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7) - mean( c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)  )

pew_it<-rbind(pew_it_a, pew_it_b)
cor(pew_it_b$X1, pew_it_b$score_a)
pew_it$country<- "Italy"


db_it_a<-data.frame(agg_zeta_it= agg_zeta_it_a$X1,agg_beta_it= agg_beta_it, names=pages_names_it )
db_it_a$State<- "State H"
db_it_b<-data.frame(agg_zeta_it= agg_zeta_it_b$X1,agg_beta_it= agg_beta_it, names=pages_names_it )
db_it_b$State<- "State L"

db_it<- rbind(db_it_a,db_it_b)
db_it$State<-factor(db_it$State,  levels = c("State L", "State H"))


p2<-ggplot(db_it, aes(x = agg_zeta_it, y = agg_beta_it, label = names)) +facet_wrap(~State, scales = "free_y")
p2<- p2+ geom_point(col = "#ff420f", size = 1)
p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
#p2<- p2 + geom_text_repel(data = db_it, size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Italy", col = "white")
p2<- p2 + geom_text_repel(data = db_it[db_it$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", title = "Italy", col = "white")
p2 <- p2 + theme(panel.grid = element_blank())  
p2 <- p2 + theme_minimal() + theme(strip.placement = "outside",legend.position="none", plot.title = element_text(face = "italic",
                                                                                                                 size = 14),
                                   strip.text.x = element_blank(),
                                   
                                   axis.title = element_text( face = "italic",size = rel(1)),
                                   axis.title.y = element_text(size = 14, angle=90,vjust =2),
                                   axis.title.x = element_blank(),
                                   axis.text=element_text(size=12),
                                   axis.line = element_blank(),
                                   axis.ticks = element_line(),
                                   panel.grid = element_blank(),
                                   panel.border = element_rect(colour = "black", fill = NA))






p2_it<-p2



result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)


gat1<-gather(result3[result3$ite %in% seq(25000, 35000, 10 ),1:4], key = "parameter", value = "value")

result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)


gat2<-gather(result8[result8$ite %in% seq(25000, 35000, 1 ),1:2], key = "parameter", value = "value")
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


gat_it<- gat
gat_it$country <- "Italy"

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




w_it<-z 



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
xi_it_agg_area$lead[length(xi_it_agg_area$lead)]<- xi_it_agg$t[length(xi_it_agg$t)] +1

q<- ggplot(xi_it_agg)+labs(x ="time", y = "State", colour = "State")
#q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
q<- q + geom_rect(data = xi_it_agg_area , aes(xmin = t, xmax = lead,  ymax  = state1),  ymin= min(as.numeric(xi_it_agg$state1)), fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
#q<- q + geom_step(aes(x = t, y = as.numeric(State)), color ="black", shape = 15 , size = 0.2)

# q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
#                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# 
q<- q + labs(title = "Italy",
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



q_it<- q


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
  geom_tile(color = "white")+ labs(x= "States", y ="States", title = "Italy" )+
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



ggheatmap_it <- ggheatmap1


##############SPAIN#####################

load("Data/Dynamic/DataEnv_SP_all.RData")
load("Code/06-Dynamic/Results/RESULT_SC_SPz_fe.RData")

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


pages_names_sp<- c(unique(EL_x$i)) # a vector of page names
# pages_names_sp<- pages_names_sp[!pages_names_sp %in% "La Voz de Asturias"]

top_names<- c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" )

beta_sp<- result[[1]]

zeta_sp_a<- result[[3]]
zeta_sp_a<-data.frame(zeta_sp_a)

zeta_sp_b<- result[[4]]
zeta_sp_b<-data.frame(zeta_sp_b)


colnames(zeta_sp_a) <-c("X1","X2", "i", "sp")
colnames(zeta_sp_b) <-c("X1","X2", "i", "sp")


zeta_sp_a_sub <-zeta_sp_a[zeta_sp_a$sp %in% seq(20000, 35000, 1), ]
zeta_sp_b_sub <-zeta_sp_b[zeta_sp_b$sp %in% seq(20000, 35000, 1), ]


zi_a_it<-data.frame(zeta_sp_a)
colnames(zi_a_it)<-c("value", "value2", "i", "it")
zi_a_it<-zi_a_it[,-2]

zi_a_it<-reshape2::dcast(zi_a_it, it ~ i)
zi_a_it<-zi_a_it[,-1]

zi_b_it<-data.frame(zeta_sp_b)
colnames(zi_b_it)<-c("value", "value2", "i", "it")
zi_b_it<-zi_b_it[,-2]

zi_b_it<-reshape2::dcast(zi_b_it, it ~ i)
zi_b_it<-zi_b_it[,-1]


# obj1<-update_zi(as.matrix(zi_a_it), as.matrix(zi_b_it))
# zi_a_it<-obj1[[1]]
# zi_b_it<-obj1[[2]]

agg_zeta_sp_a<-aggregate(zeta_sp_a_sub[,1:2], by = list(zeta_sp_a_sub$i), FUN = mean )
agg_zeta_sp_a$X1<-colMeans(zi_a_it[seq(20000, 35000, ), ])

agg_beta_sp<-colMeans(beta_sp[seq(20000, 35000, 1), ])
agg_zeta_sp_a$names<- pages_names_sp

agg_zeta_sp_b<-aggregate(zeta_sp_b_sub[,1:2], by = list(zeta_sp_b_sub$i), FUN = mean )
agg_zeta_sp_b$X1<-colMeans(zi_b_it[seq(20000, 35000, ), ])

agg_beta_sp<-colMeans(beta_sp[seq(20000, 35000, 1), ])
agg_zeta_sp_b$names<- pages_names_sp


pew_sp_a<- agg_zeta_sp_a[agg_zeta_sp_a$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" ), ]
pew_sp_a$State<- "State H"
pew_sp_a$score_a<-  c(4.5, 3.8, 4.2, 3.4, 3.5) - mean( c(4.5, 3.8, 4.2, 3.4, 3.5))
pew_sp_a$score_b<-  c(3.3, 3.1, 3.2, 3.1, 2.8) - mean( c(3.3, 3.1, 3.2, 3.1, 2.8) )


pew_sp_b<- agg_zeta_sp_b[agg_zeta_sp_b$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" ), ]
pew_sp_b$State<- "State L"
pew_sp_b$score_a<-  c(4.5, 3.8, 4.2, 3.4, 3.5) - mean( c(4.5, 3.8, 4.2, 3.4, 3.5))
pew_sp_b$score_b<-  c(3.3, 3.1, 3.2, 3.1, 2.8) - mean( c(3.3, 3.1, 3.2, 3.1, 2.8) )

pew_sp<-rbind(pew_sp_a, pew_sp_b)

pew_sp$country<- "Spain"


db_sp_a<-data.frame(agg_zeta_sp= agg_zeta_sp_a$X1,agg_beta_sp= agg_beta_sp, names=pages_names_sp )
db_sp_a$State<- "State H"
db_sp_b<-data.frame(agg_zeta_sp= agg_zeta_sp_b$X1,agg_beta_sp= agg_beta_sp, names=pages_names_sp )
db_sp_b$State<- "State L"

db_sp<- rbind(db_sp_a,db_sp_b)
db_sp$State<-factor(db_sp$State,  levels = c("State L", "State H"))

p2<-ggplot(db_sp, aes(x = agg_zeta_sp, y = agg_beta_sp, label = names)) +facet_wrap(~State, scales = "free_y")
p2<- p2+ geom_point(col = "#ff420f", size = 1)
p2<- p2+ scale_x_continuous(limits = symmetric_limits) +  scale_y_continuous(limits = symmetric_limits) 
#p2<- p2 + geom_text_repel(data = db_sp,size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y ="Individual Effect", title = "Spain", col = "whspe")
p2<- p2 + geom_text_repel(data = db_sp[db_sp$names %in% top_names,],size = 3,  min.segment.length = 0)  + labs( x = "Latent Leaning", y ="Individual Effect", title = "Spain", col = "whspe")
p2 <- p2 + theme(panel.grid = element_blank())  
p2 <- p2 + theme_minimal() + theme(strip.placement = "outside", legend.position="none", plot.title = element_text(face = "italic",
                                                                                                                  size = 14),
                                   strip.text.x = element_blank(),
                                   
                                   axis.title = element_text( face = "italic",size = rel(1)),
                                   axis.title.y = element_text(size = 14, angle=90,vjust =2),
                                   axis.title.x = element_text(size = 14, vjust = -0.2),
                                   axis.text=element_text(size=12),
                                   axis.line = element_blank(),
                                   axis.ticks = element_line(),
                                   panel.grid = element_blank(),
                                   panel.border = element_rect(colour = "black", fill = NA))



p2_sp<-p2




result3<- result[[2]]

colnames(result3)<-c("phi", "gamma_0", "gamma_1", "tau", "ite")
result3<-data.frame(result3)


gat1<-gather(result3[result3$ite %in% seq(25000, 35000, 10 ),1:4], key = "parameter", value = "value")

result8<-data.frame(result$Sigma_z_ite)
result8$X1<- result8$X1^2
result8$X2<- result8$X2^2
result8$sa_mean <- cumsum(result8$X1)/(1:length(result8$X1))
result8$sb_mean <- cumsum(result8$X2)/(1:length(result8$X2))
result8$ite<- 1:length(result8$X2)


gat2<-gather(result8[result8$ite %in% seq(25000, 35000, 1 ),1:2], key = "parameter", value = "value")
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

gat_sp<- gat
gat_sp$country <- "Spain"


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
z <- z + labs(title = "Spain", 
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




w_sp<-z 



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
xi_it_agg_area$lead[length(xi_it_agg_area$lead)]<- xi_it_agg$t[length(xi_it_agg$t)] +1

q<- ggplot(xi_it_agg)+labs(x ="time", y = "State", colour = "State")
#q<- q + geom_point(aes(x = t, y = State), color ="#ff420f", shape = 15 , size = 0.1)
q<- q + geom_rect(data = xi_it_agg_area , aes(xmin = t, xmax = lead,  ymax  = state1),  ymin= min(as.numeric(xi_it_agg$state1)), fill = "#ff420f",  alpha = 0.1, inherit.aes = F)
q<- q + geom_point(data = xi_it_agg , aes(x = t,  y= as.numeric(State)), size = 0.6, col = "#ff420f", alpha = 0.3, inherit.aes = F)
#q<- q + geom_step(aes(x = t, y = as.numeric(State)), color ="black", shape = 15 , size = 0.2)

# q<- q + geom_rect(aes(xmin = t, xmax = lead(t), 
#                       ymin = 0, ymax  = 1-(as.numeric(State)-1)), col ="red", alpha = 0.01)
# 
q<- q + labs(title = "Spain",
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



q_sp<- q



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
  geom_tile(color = "white")+ labs(x= "States", y ="States", title = "Spain" )+
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



ggheatmap_sp <- ggheatmap1


####
PEW_all<-rbind(pew_de, pew_fr, pew_it, pew_sp)
PEW_all$State<- factor(PEW_all$State, levels = c("State L", "State H")  )

##########

w<- (w_fr|w_de)/(w_it|w_sp) & theme(axis.text=element_text(size= 6), axis.text.x = element_text(angle = 30,vjust = 0.7))
w

ggsave(w, filename = "Figures/Dynamic/Figure8.pdf", units = "cm", width = 16*2, height = 9*2 )


##########
#figure 8 v2

gat_all<-rbind(gat_fr, gat_de, gat_it, gat_sp)

sturges <- function(x){ pretty(range(x),
                               n = nclass.Sturges(x),
                               min.n = 1)}

z <- ggplot(gat_all,  aes(x=value, group = country, linetype = country )) +  geom_density(fill = "#ff420f", alpha = 0.3)  + facet_wrap(.~level, scales = "free",  labeller = label_parsed)
z<- z + geom_hline(yintercept = 0 , col = "white", linewidth = 1)
# z <-  z+ geom_histogram(data = gat[gat$level == "gamma[0]",], aes( y = after_stat(density)), binwidth=0.1, fill = "#ff420f", alpha = 0.6 , col ="black", breaks = sturges(gat[gat$level == "gamma[0]","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "gamma[1]",],  aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "gamma[1]","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "phi",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "phi","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "sigma[L]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[L]^2","value"] ) ) 
# z <- z+ geom_histogram(data = gat[gat$level == "sigma[H]^2",],   aes( y = after_stat(density)), binwidth=0.01, fill = "#ff420f", alpha = 0.6 , col ="black" ,  breaks = sturges(gat[gat$level == "sigma[H]^2","value"] ) ) 
# z<- z + geom_density(linetype = 1,  n = 10000, adjust = 2)
z <- z + scale_linetype_manual(values=c("solid", "twodash", "dotted", "dashed")) 
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "phi") , fun = dgamma, args = list(shape=0.01, rate =0.01), geom = "area", fill ="black", alpha = 0.3, linetype = "dashed", colour = "red", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[0]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3, linetype = "dashed", colour = "red", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "gamma[1]"),fun = dnorm, args = list(mean = 0, sd = 15^2), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[L]^2"),fun = dinvgamma, args = list(shape = 0.01 ,scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)
z<-z+  stat_function( data = data.frame(value = 0, yaxis = 0, level = "sigma[H]^2"),fun = dinvgamma, args = list(shape = 0.01, scale=0.01), geom = "area",fill ="black", alpha = 0.3,linetype = "dashed", colour = "red", inherit.aes = F)
#z<-z+ facet_wrap( .~ factor(level ,  levels = c("gamma[0]", "gamma[1]", "phi", "sigma[L]^2", "sigma[H]^2")), ncol = 6, scales = "free",labeller = label_parsed)
# z<- z + scale_shape_manual(labels =  parse_format())
z <- z + labs(x = "Value",  y = "", linetype="Country")+  theme_minimal() + theme( legend.direction = "vertical", plot.title = element_text(face = "italic",
                                                                                                                size = 22, hjust = 0.5),
                                                              legend.key.size = unit(1.5, 'cm'),
                                                              legend.title = element_text(size=20), #change legend title font size
                                                              legend.text = element_text(size=18),
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




w_all<-z + guides(linetype=guide_legend(ncol=2))


library(gtable)
library(cowplot)

shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

library(grid)


## Initiate writing to PDF file
pdf("Figures/Dynamic/Figure8_v2.pdf",  width = 16, height = 9)

## Create a graphical object g here
grid.draw(shift_legend(w_all))
# print it

## Stop writing to the PDF file
dev.off()

##########


pp<- (p2_fr/p2_de/p2_it/p2_sp)  & theme(axis.title = element_text(size = 16) )
pp<- (wrap_ggplot_grob(gp)/p2_de/p2_it/p2_sp)  & theme(axis.title = element_text(size = 16) )


ggsave(pp, filename = "Figures/Dynamic/Figure9.pdf", units = "cm", width = 21, height = 29.7 )


p5<- ggplot(PEW_all) + facet_wrap(.~State) +geom_point(aes(x =  X1 , y = score_a ), col = "#ff420f") +geom_text_repel(aes(x =  X1 , y = score_a, label = names ), size = 3) + theme_minimal()

p5<- p5 + labs(x = "Latent Leaning", y = "PEW score") # + geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
p5<- p5 + stat_cor(aes(x =  X1 , y = score_a, label = after_stat(r.label))) #+ xlim(-0.25, 0.25)  
p5 <- p5 + theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                       size = 16),
                                   strip.text.x = element_text(size = 16,face = "italic"),
                                   
                                   axis.title = element_text( face = "italic",size = rel(1)),
                                   axis.title.y = element_text(size = 16, angle=90,vjust =2),
                                   axis.title.x = element_text(size = 16, vjust = -0.2),
                                   axis.text=element_text(size=14),
                                   axis.line = element_blank(),
                                   axis.ticks = element_line(),
                                   panel.grid = element_blank(),
                                   legend.title = element_text(face="italic", size=16),
                                   legend.text = element_text(size=14),
                                   panel.border = element_rect(colour = "black", fill = NA))





p5


ggsave(p5, filename = "Figures/Dynamic/Figure10.pdf", units = "cm", width = 16*2, height = 9*2 )



q_all <- (q_fr|q_de)/(q_it|q_sp)+(ggheatmap_fr|ggheatmap_de|ggheatmap_it|ggheatmap_sp) + plot_annotation(tag_levels = 'A') & theme(plot.title = element_text(size = 24), axis.text = element_text(size = 18)  )
ggsave(q_all, filename = "Figures/Dynamic/Figure11.pdf", units = "cm", width = 16*3, height = 9*3 )


# ####### DIC ############
# #load the libraries
# 
# library(LaplacesDemon)
# library(mvtnorm)
# library(reshape2)
# library(plyr)
# library(tidyr)
# library(data.table)
# library(MASS)
# library(vegan)
# library(Rcpp)
# library(RcppArmadillo)
# 
# setwd("~/Desktop/political_programs/political_programs/01_IT")
# sourceCpp("Rcpp_rf_dynamic.cpp")
# 
# logistic_f<-function(x){(1/(1+exp(-1*x)))} #logistic function
# 
# 
# M<- as.matrix(result[[7]])
# Time <- 731
# 
# result1<- result[[1]]
# beta = colMeans(result1[seq(25000, 35000, 1), ])
# phi <- mean(result3$phi[seq(25000, 35000, 1)])
# gamma0 <- mean(result3$gamma_0[seq(25000, 35000, 1)])
# gamma1 <- mean(result3$gamma_1[seq(25000, 35000, 1)])
# 
# beta_v = vec_nbet_tri_a( beta,  Npages)
# beta_v_not = vec_nbet_tri( beta,  Npages)
# 
# result4<- result[[3]]
# result5<- result[[4]]
# 
# colnames(result4) <-c("X1","X2", "i", "it")
# result4<-data.frame(result4)
# 
# colnames(result5) <-c("X1","X2", "i", "it")
# result5<-data.frame(result5)
# 
# 
# zi_b_it_sub<-result5[result5$it %in% seq(25000, 35000, 1), ]
# agg_zi_b_it<-aggregate(zi_b_it_sub[,1:2], by = list(zi_b_it_sub$i), FUN = mean)
# 
# 
# zi_a_it_sub<-result4[result4$it %in%seq(25000, 35000, 1), ]
# agg_zi_a_it<-aggregate(zi_a_it_sub[,1:2], by = list(zi_a_it_sub$i), FUN = mean )
# 
# 
# #get the distance between pages in state a
# A<-matrix(0, Npages,Npages)
# ui<-which(upper.tri(A, diag = F), arr.ind=T)
# upper_index<- paste0(ui[,1],"_",ui[,2])
# 
# eucl_a<- as.matrix(dist(agg_zi_a_it[,2], method = "euclidean", diag = T, upper = F, p = 2))
# eucl_a <- reshape2::melt(eucl_a)
# colnames(eucl_a) <- c("j", "i", "d")
# eucl_a$i_j<- paste0(eucl_a$i,"_",eucl_a$j)
# eucl_a <- eucl_a[eucl_a$i_j %in%  upper_index, ]
# eucl_a<-eucl_a[order(eucl_a$i, eucl_a$j), ]
# 
# #get the distance between pages in state b
# 
# eucl_b<- as.matrix(dist(agg_zi_b_it[,2], method = "euclidean", diag = T,  p = 2))
# eucl_b <- reshape2::melt(eucl_b)
# colnames(eucl_b) <- c("j", "i", "d")
# eucl_b$i_j<- paste0(eucl_b$i,"_",eucl_b$j)
# eucl_b <- eucl_b[eucl_b$i_j %in%  upper_index, ]
# eucl_b<-eucl_b[order(eucl_b$i, eucl_b$j), ]
# 
# #compute likelihoods for the dynamic edge set
# 
# EL_princ$distance_a<- eucl_a$d
# EL_princ$distance_b<-eucl_b$d
# 
# EL_princ$lambda_a<- exp( beta_v  + beta_v_not -1*(EL_princ$distance_a)  )
# EL_princ$lambda_b <- exp(beta_v  + beta_v_not -1*(EL_princ$distance_b) )
# 
# EL_princ$g_a<- dpois(EL_princ$w, lambda = EL_princ$lambda_a, log = T)  
# EL_princ$g_b<- dpois(EL_princ$w, lambda = EL_princ$lambda_b, log = T)
# 
# data_t = as.data.table(EL_princ[,c("t","g_a", "g_b")] )
# 
# ProdG =  data_t[,list( prodG_a = sum(g_a), prodG_b = sum(g_b)),  by = "t"]
# 
# #compute likelihoods for the dynamic plane
# 
# DBplane$X1a <- agg_zi_a_it$X1
# DBplane$X1b <-  agg_zi_b_it$X1
# 
# DBplane$a_a<- logistic_f(gamma0 + gamma1*DBplane$X1a)*phi
# DBplane$b_a<- (1- logistic_f( gamma0 + gamma1*DBplane$X1a))*phi
# 
# DBplane$a_b<- logistic_f(gamma0 + gamma1*DBplane$X1b)*phi
# DBplane$b_b<- (1- logistic_f( gamma0 + gamma1*DBplane$X1b))*phi
# 
# DBplane$B_a <- dbeta(DBplane$leaning, shape1 = DBplane$a_a, shape2 =  DBplane$b_a, log =T)
# DBplane$B_b <- dbeta(DBplane$leaning, shape1 = DBplane$a_b, shape2 =  DBplane$b_b,  log =T)
# 
# data_t = as.data.table(DBplane[,c("t","B_a", "B_b")] )
# ProdH = data_t[,list( prodB_a = sum(B_a), prodB_b = sum(B_b) ), by = 't']
# 
# #compute ETA as in Hamilton but in log scale
# 
# Eta<-data.frame(t = 1:Time, likelihood_a = ProdG$prodG_a + ProdH$prodB_a ,
#                 likelihood_b = ProdG$prodG_b + ProdH$prodB_b )
# 
# result7<- result[[6]]
# colnames(result7)<-c("state1", "state2", "t", "ite")
# result7<-data.frame(result7)
# 
# library(ggplot2)
# library(lubridate)
# xi_it_sub<-result7[result7$ite %in% seq(25000, 35000, 1) ,]
# xi_it_agg<-aggregate(xi_it_sub[,1:2], by = list(t = xi_it_sub$t), FUN = median)
# 
# #compute loglik bar
# 
# likeBar<- sum(Eta$likelihood_a*xi_it_agg$state1 + Eta$likelihood_b*(1-xi_it_agg$state1) )
# 
# result8<- as.matrix(result[[7]])
# 
# DIC<- -2*mean(result8[seq(25000, 35000, 1), 3]) +2*(likeBar - mean(result8[seq(25000, 35000, 1), 3]) )
# DIC
# 
# 
# 
# ####### TRACE PLOTS ############
# 
# country <- "sp"
# 
# result1<-result[[1]]
# result1<-data.frame(result1)
# result1$ite<-1:nrow(result1)
# result1$X1_mean <- cumsum(result1$X1)/(1:length(result1$X1))
# result1$X11_mean <- cumsum(result1$X11)/(1:length(result1$X11))
# 
# 
# 
# p10<-ggplot(result1) + geom_line(aes(x = ite , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = X1_mean), col = "black") 
# p10<- p10 +  labs(title = expression(paste("Parameter ", alpha[1])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                                size = rel(1.2), hjust = 0.5),
#                                                                                                                            
#                                                                                                                            axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                            axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                            axis.title.x = element_text(vjust = -0.2),
#                                                                                                                            axis.line = element_line(colour="black"),
#                                                                                                                            axis.ticks = element_line(),
#                                                                                                                            legend.title = element_text(face="italic"),
#                                                                                                                            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p10
# 
# 
# p11<-ggplot(result1) + geom_line(aes(x = ite , y = X11), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = X11_mean), col = "black") 
# p11<- p11 +  labs(title = expression(paste("Parameter ", alpha[11])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                                 size = rel(1.2), hjust = 0.5),
#                                                                                                                             
#                                                                                                                             axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                             axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                             axis.title.x = element_text(vjust = -0.2),
#                                                                                                                             axis.line = element_line(colour="black"),
#                                                                                                                             axis.ticks = element_line(),
#                                                                                                                             legend.title = element_text(face="italic"),
#                                                                                                                             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p11
# 
# 
# result3$phi_mean <- cumsum(result3$phi)/(1:length(result3$phi))
# result3$gamma0_mean <- cumsum(result3$gamma_0)/(1:length(result3$gamma_0))
# result3$gamma1_mean <- cumsum(result3$gamma_1)/(1:length(result3$gamma_1))
# 
# 
# p7<-ggplot(result3) + geom_line(aes(x = ite , y = phi), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = phi_mean), col = "black")  
# p7<- p7 +  labs(title = expression(paste("Parameter ", phi)), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                         size = rel(1.2), hjust = 0.5),
#                                                                                                                     
#                                                                                                                     axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                     axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                     axis.title.x = element_text(vjust = -0.2),
#                                                                                                                     axis.line = element_line(colour="black"),
#                                                                                                                     axis.ticks = element_line(),
#                                                                                                                     legend.title = element_text(face="italic"),
#                                                                                                                     strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p7
# 
# 
# p8<-ggplot(result3) + geom_line(aes(x = ite , y = gamma_0), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = gamma0_mean), col = "black")  
# p8<- p8 +  labs(title = expression(paste("Parameter ", gamma[0])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                              size = rel(1.2), hjust = 0.5),
#                                                                                                                          
#                                                                                                                          axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                          axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                          axis.title.x = element_text(vjust = -0.2),
#                                                                                                                          axis.line = element_line(colour="black"),
#                                                                                                                          axis.ticks = element_line(),
#                                                                                                                          legend.title = element_text(face="italic"),
#                                                                                                                          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p8
# 
# 
# p9<-ggplot(result3) + geom_line(aes(x = ite , y =  gamma_1 ), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = ite , y = gamma1_mean), col = "black")
# p9<- p9 +  labs(title = expression(paste("Parameter ", gamma[1])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                              size = rel(1.2), hjust = 0.5),
#                                                                                                                          
#                                                                                                                          axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                          axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                          axis.title.x = element_text(vjust = -0.2),
#                                                                                                                          axis.line = element_line(colour="black"),
#                                                                                                                          axis.ticks = element_line(),
#                                                                                                                          legend.title = element_text(face="italic"),
#                                                                                                                          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p9
# 
# 
# mc2<-(p10|p11)/(p8|p9)/(p7)  + plot_annotation(tag_levels = 'A')
# ggsave(mc2, filename = paste("trace_plane_", country, ".pdf"), units = "cm", width = 16*2, height = 9*2 )
# 
# 
# 
# library(patchwork)
# result4<-result[[3]]
# result4<-data.frame(result4)
# colnames(result4) <-c("X1","X2", "i", "it")
# try<-data.frame(result4)
# 
# 
# 
# try = try[try$i == 20, ]
# hist(try$X1[20000:35000])
# plot(try$X1[20000:35000], type = "l")
# 
# 
# try$mean_z1 <- cumsum(try$X31)/(1:length(try$X31))
# 
# 
# p12<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  
# p12<- p12 +  labs(title = expression(paste("Parameter ", zeta["1, H"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                                    size = rel(1.2), hjust = 0.5),
#                                                                                                                                
#                                                                                                                                axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                                axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                                axis.title.x = element_text(vjust = -0.2),
#                                                                                                                                axis.line = element_line(colour="black"),
#                                                                                                                                axis.ticks = element_line(),
#                                                                                                                                legend.title = element_text(face="italic"),
#                                                                                                                                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p12
# 
# 
# 
# try<-data.frame(result4)
# try = try[try$i == 11, ]
# 
# try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))
# 
# 
# p13<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black") 
# p13<- p13 +  labs(title = expression(paste("Parameter ", zeta["11, H"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                                     size = rel(1.2), hjust = 0.5),
#                                                                                                                                 
#                                                                                                                                 axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                                 axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                                 axis.title.x = element_text(vjust = -0.2),
#                                                                                                                                 axis.line = element_line(colour="black"),
#                                                                                                                                 axis.ticks = element_line(),
#                                                                                                                                 legend.title = element_text(face="italic"),
#                                                                                                                                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p13
# 
# 
# 
# 
# 
# 
# 
# 
# result5<-result[[4]]
# result5<-data.frame(result5)
# colnames(result5) <-c("X1","X2", "i", "it")
# try<-data.frame(result5)
# try = try[try$i == 20, ]
# 
# 
# 
# try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))
# 
# 
# p14<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  
# p14<- p14 +  labs(title = expression(paste("Parameter ", zeta["1, L"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                                    size = rel(1.2), hjust = 0.5),
#                                                                                                                                
#                                                                                                                                axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                                axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                                axis.title.x = element_text(vjust = -0.2),
#                                                                                                                                axis.line = element_line(colour="black"),
#                                                                                                                                axis.ticks = element_line(),
#                                                                                                                                legend.title = element_text(face="italic"),
#                                                                                                                                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p14
# 
# 
# 
# try<-data.frame(result5)
# try = try[try$i == 11, ]
# try$mean_z1 <- cumsum(try$X1)/(1:length(try$X1))
# 
# p15<-ggplot(try) + geom_line(aes(x = it , y = X1), col = "#ff8566", alpha = 0.8)+geom_line(aes(x = it , y = mean_z1), col = "black")  
# p15<- p15 +  labs(title = expression(paste("Parameter ", zeta["11, L"])), x = "Iteration", y = "Value")+ theme_classic()+ theme(legend.position="bottom", plot.title = element_text(face = "italic",
#                                                                                                                                                                                     size = rel(1.2), hjust = 0.5),
#                                                                                                                                 
#                                                                                                                                 axis.title = element_text(face = "italic",size = rel(1)),
#                                                                                                                                 axis.title.y = element_text(angle=90,vjust =2),
#                                                                                                                                 axis.title.x = element_text(vjust = -0.2),
#                                                                                                                                 axis.line = element_line(colour="black"),
#                                                                                                                                 axis.ticks = element_line(),
#                                                                                                                                 legend.title = element_text(face="italic"),
#                                                                                                                                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0") )
# p15
# 
# 
# 
# 
# library(patchwork)
# mc3<-(p12|p13)/(p14|p15)  + plot_annotation(tag_levels = 'A')
# ggsave(mc3, filename = paste("trace_zeta_", country, ".pdf"), units = "cm", width = 16*2, height = 9*2 )
# 
# 
# 
