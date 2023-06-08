

rm(list = ls())

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
# check<-reshape(zeta_de, idvar = "it", timevar = "i", direction = "wide")
# check<-check[,-1]
# plot(c(as.matrix(check[2000:15000,])), c(beta_de[2000:15000,]), pch = ".")


top_names<-  c("Bild", "FAZ.NET - Frankfurter Allgemeine Zeitung", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" )

pew_de<- agg_zeta_de[agg_zeta_de$names %in% c("Bild", "FAZ.NET - Frankfurter Allgemeine Zeitung", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" ), ]
pew_de$score_a<- c(3.6, 3.2, 3.2, 2.7, 3)
pew_de$score_b<- c(3.1, 2.9, 3, 2.8, 2.8)

mod<-lm(pew_de$score_a ~ pew_de$X1)
plot(pew_de$score_a ~ pew_de$X1)
abline(mod$coefficients[1],mod$coefficients[2])
plot(log(DBplane$ncomments), colMeans(beta_de[10000:15000,]))
cor(log(DBplane$ncomments), colMeans(beta_de[10000:15000,]))
lc<-log(DBplane$ncomments)
mod2<-lm(  colMeans(beta_de[10000:15000,]) ~ lc)
abline(mod2$coefficients[1], mod2$coefficients[2])


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
# gat_prior<-gat
# gat_prior$value <-c(prior_phi, prior_g0, prior_g1)
# gat_prior$type <- "prior"
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
#pew_de$X1<- pew_de$X1 -  mean(pew_de$X1)
pew_de$score_a<- pew_de$score_a -  mean(pew_de$score_a)


#


result4<- result[[3]]

result4<-data.frame(result4)
colnames(result4) <-c("X1","X2", "i", "it")
try<-data.frame(result4)


result4_t = data.table(result4)
distance_a<- result4_t[ !i %in% c(15, 36) ,mean(as.numeric(dist(X1))), by = it]

gat3<-data.frame(dh = distance_a$V1)
colnames(gat3)<-c("value")


p <- ggplot(gat3,  aes(x = 0 , y=value)) + stat_summary(fun.data= mean_sdl, fun.args = list(mult=3), 
                                                        geom="pointrange", color="#ff420f")
p<- p+ scale_shape_manual(labels =  parse_format())
p<- p + labs( title = "Germany",y = "Avg. Distance")+  theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                                                                 size = 36, hjust = 0.5),
                                                                               axis.text=element_text(size=16),
                                                                               strip.text.x = element_text(size = 18, face = "italic"),
                                                                               axis.title.x = element_blank(),
                                                                               axis.text.x = element_blank(),
                                                                               axis.title.y = element_text(size = 22, angle=90,vjust =2, face = "italic"),
                                                                               axis.line = element_blank(),
                                                                               axis.ticks = element_line(),
                                                                               axis.ticks.x = element_blank(),
                                                                               panel.grid = element_blank(),
                                                                               panel.border = element_rect(colour = "black", fill = NA))



p_de <- p


############

load("Code/05-Static/RESULT_single_fr.RData")
load("Data/Static/Data_Env_single_FR.RData")
pages_names_fr<- c(unique(EL_x$i)) # a vector of page names

top_names<-c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart")

beta_fr<- result[[1]]

zeta_fr<- result[[3]]
zeta_fr<-data.frame(zeta_fr)
colnames(zeta_fr) <-c("X1","X2", "i", "it")

zeta_fr_sub <-zeta_fr[zeta_fr$it %in% seq(2000, 15000, 1), ]

agg_zeta_fr<-aggregate(zeta_fr_sub[,1:2], by = list(zeta_fr_sub$i), FUN = mean )
agg_beta_fr<-colMeans(beta_fr[seq(2000, 15000, 1), ])
agg_zeta_fr$names<- pages_names_fr


zeta_fr_sub <-zeta_fr[zeta_fr$it %in% seq(2000, 15000, 1), ]
zeta_fr<-zeta_fr[,-2]
# check<-reshape(zeta_fr, idvar = "it", timevar = "i", direction = "wide")
# check<-check[,-1]
# plot(c(as.matrix(check[2000:15000,])), c(beta_fr[2000:15000,]), pch = ".")
# 

db_fr<-data.frame(agg_zeta_fr= agg_zeta_fr$X1,agg_beta_fr= agg_beta_fr, names=pages_names_fr )
p2<-ggplot(db_fr, aes(x = agg_zeta_fr, y = agg_beta_fr, label = names))
p2<- p2+ geom_point(col = "#ff420f", size = 1) 
p2<- p2 + geom_text_repel(data = db_fr[db_fr$names %in% top_names,] ,size = 2,  min.segment.length = 0)  + labs(  x = "Latent Leaning", y = "Individual Effect", title = "France", col = "white")
p2 <- p2 + theme(panel.grid = element_blank())  + xlim(-2.5,2.5) +ylim(-2,7)
p2 <- p2 + theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
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



p2


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




w_fr<-z +  coord_cartesian_panels(
  panel_limits = tibble::tribble(
    ~level, ~xmin, ~xmax,
    "gamma[0]"     ,   mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]),
    "gamma[1]"     , mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]),
    "phi"     ,   0,     mean(gat$value[gat$parameter == "phi" & gat$type == "posterior" ]) + 7*sd(gat$value[gat$parameter == "phi" & gat$type == "posterior" ])
    
  ))


pew_fr<- agg_zeta_fr[agg_zeta_fr$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
pew_fr$score_a<- c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1)
pew_fr$score_b<- c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3)


#pew_fr$X1<- pew_fr$X1 -  mean(pew_fr$X1)
pew_fr$score_a<- pew_fr$score_a -  mean(pew_fr$score_a)
#########


load("Code/05-Static/RESULT_single_it.RData")
load("Data/Static/Data_Env_single_IT.RData")
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



pew_it<- agg_zeta_it[agg_zeta_it$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ),]
pew_it$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5)
pew_it$score_b<-c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)

#pew_it$X1<- pew_it$X1 -  mean(pew_it$X1)
pew_it$score_a<- pew_it$score_a -  mean(pew_it$score_a)

##########



load("Code/05-Static/RESULT_single_sp.RData")
load("Data/Static/Data_Env_single_SP.RData")
pages_names_sp<- c(unique(EL_x$i)) # a vector of page names

beta_sp<- result[[1]]

zeta_sp<- result[[3]]
zeta_sp<-data.frame(zeta_sp)
colnames(zeta_sp) <-c("X1","X2", "i", "it")

zeta_sp_sub <-zeta_sp[zeta_sp$it %in% seq(2000, 15000, 1), ]

agg_zeta_sp<-aggregate(zeta_sp_sub[,1:2], by = list(zeta_sp_sub$i), FUN = mean )
agg_beta_sp<-colMeans(beta_sp[seq(2000, 15000, 1), ])
agg_zeta_sp$names<- pages_names_sp

zeta_sp_sub <-zeta_sp[zeta_sp$it %in% seq(2000, 15000, 1), ]
zeta_sp<-zeta_sp[,-2]
# check<-reshape(zeta_sp, idvar = "it", timevar = "i", direction = "wide")
# check<-check[,-1]
# plot(c(as.matrix(check[2000:15000,])), c(beta_sp[2000:15000,]), pch = ".")


top_names<- c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" )


db_sp<-data.frame(agg_zeta_sp= agg_zeta_sp$X1,agg_beta_sp= agg_beta_sp, names=pages_names_sp )
p4<-ggplot(db_sp, aes(x = agg_zeta_sp, y = agg_beta_sp, label = names))
p4<- p4+ geom_point(col = "#ff420f", size = 1) 
p4<- p4 + geom_text_repel(data = db_sp[db_sp$names %in% top_names,] ,size = 2,  min.segment.length = 0)  + labs( x = "Latent Leaning", y = "Individual Effect", y = "FE", title = "Spain", col = "white")
p4 <- p4 + theme(panel.grid = element_blank())  + xlim(-2.5,2.5) +ylim(-2,7)
p4 <- p4 + theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
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



p4




pew_sp<- agg_zeta_sp[agg_zeta_sp$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" ),]
pew_sp$X1<- pew_sp$X1 -  mean(pew_sp$X1)
pew_sp$score_a<- c(4.5, 3.8, 4.2, 3.4, 3.5)
pew_sp$score_b<- c(3.3, 3.1, 3.2, 3.1, 2.8)




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




w_sp<-z +  coord_cartesian_panels(
  panel_limits = tibble::tribble(
    ~level, ~xmin, ~xmax,
    "gamma[0]"     ,   mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[0]" & gat$type == "posterior" ]),
    "gamma[1]"     , mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) - 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) ,      mean(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]) + 5*sd(gat$value[gat$parameter == "gamma[1]" & gat$type == "posterior" ]),
    "phi"     ,   0,     mean(gat$value[gat$parameter == "phi" & gat$type == "posterior" ]) + 7*sd(gat$value[gat$parameter == "phi" & gat$type == "posterior" ])
    
  ))


#pew_sp$X1<- pew_sp$X1 -  mean(pew_sp$X1)
pew_sp$score_a<- pew_sp$score_a -  mean(pew_sp$score_a)

PEW_all<-rbind(pew_de, pew_fr, pew_it, pew_sp)

##########


p5<- ggplot(PEW_all) +geom_point(aes(x =  X1 , y = score_a ), col = "#ff420f") +geom_text_repel(aes(x =  X1 , y = score_a, label = names ), size = 3) + theme_minimal()
p5<- p5 + labs(x = "Latent Leaning", y = "PEW score") # + geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
p5<- p5 + stat_cor(aes(x =  X1 , y = score_a, label = ..r.label..)) #+ xlim(-0.25, 0.25)  
p5 <- p5 + theme_minimal() + theme(legend.position="bottom", plot.title = element_text(face = "italic",
                                                                                       size = 22, hjust = 0.5),
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





p5


p2<-p2 + theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title = element_text(face = "italic",size = 16, hjust = 0.5),axis.text=element_text(size=10))
p1<-p1 + theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title = element_text(face = "italic",size = 16, hjust = 0.5),axis.text=element_text(size=10))
p3<-p3 + theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title = element_text(face = "italic",size = 16, hjust = 0.5),axis.text=element_text(size=10))
p4<-p4 + theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title = element_text(face = "italic",size = 16, hjust = 0.5),axis.text=element_text(size=10))
p5<- p5 + labs(title = "PEW Score - Correlation") + theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), plot.title = element_text(face = "italic",size = 16, hjust = 0.5))


design<-"1255
         3455"


p<- p2 + p1 + p3 + p4 + p5+ plot_annotation(tag_levels = 'A') + plot_layout(design = design)
p

ggsave(p, filename = "Figures/Static/Figure7.pdf", units = "cm", width = 16*2, height = 9*2 )



w<- (w_fr|w_de)/(w_it|w_sp)
w

ggsave(w, filename = "Figures/Static/FigureI1.pdf", units = "cm", width = 16*2, height = 9*2 )

##########