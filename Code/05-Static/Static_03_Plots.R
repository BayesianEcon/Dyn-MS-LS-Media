 
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
library(igraph)
library(sbm)


########CHANGE YOUR PATH ###########
setwd("~/Documents/GitHub/Dyn-MS-LS-Media/")
#####################################

load("Code/05-Static/RESULT_single_de.RData")
load("Data/Static/Data_Env_single_DE.RData")
source("Code/05-Static/Misc/CoordCartesian.R")

el<-EL_princ[,1:3]
colnames(el)<-c("from","to", "weight")
g<-graph_from_data_frame(el)

A<-as.matrix(as_adjacency_matrix(g, attr= "weight"))
A= A+t(A)

# ADJ<-r_to_py(A)
# nr<-r_to_py(ncol(A))
# py$ADJ<- ADJ
# py$n<- nr
# source_python("sbm_gtools.py")
#Membership<-factor(py$membership)

# source("degcor-sbm.R")
# g<-graph_from_adjacency_matrix(A,weighted = T)
# Membership_de<-factor(Membership)

# mySimpleSBM <- estimateSimpleSBM(A, 'poisson',
#                                  estimOptions = list(plot = FALSE))

g<-graph_from_adjacency_matrix(A,weighted = T)

E(g)$weight
g<-as.undirected(g)
cl<- cluster_fast_greedy(g ,weights  = E(g)$weight )
#cl<- cluster_walktrap(g ,weights  = E(g)$weight )
#cl<- cluster_louvain(g ,weights  = E(g)$weight )

Membership<-membership(cl)
Membership_de<-factor(Membership)

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
p1<- p1+ geom_point(aes(col = Membership_de, shape = Membership_de), size = 3) 
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


#pew_de$X1<- pew_de$X1 -  mean(pew_de$X1)
pew_de$score_a<- pew_de$score_a -  mean(pew_de$score_a)

pew_de$Membership<- factor(as.numeric(Membership[agg_zeta_de$names %in% c("Bild", "FAZ.NET - Frankfurter Allgemeine Zeitung", "RTL Aktuell", "Süddeutsche Zeitung", "SPIEGEL ONLINE" )]))
pew_de$Country<- "DE"

############

load("Code/05-Static/RESULT_single_fr.RData")
load("Data/Static/Data_Env_single_FR.RData")
pages_names_fr<- c(unique(EL_x$i)) # a vector of page names

top_names<-c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart")

el<-EL_princ[,1:3]
colnames(el)<-c("from","to", "weight")
g<-graph_from_data_frame(el)

A<-as.matrix(as_adjacency_matrix(g, attr= "weight"))
A= A+t(A)

# source("degcor-sbm.R")
# Membership_fr<-Membership

# ADJ<-r_to_py(A)
# nr<-r_to_py(ncol(A))
# py$ADJ<- ADJ
# py$n<- nr
# source_python("sbm_gtools.py")
# Membership<-factor(py$membership)
# Membership_fr<-Membership

g<-graph_from_adjacency_matrix(A,weighted = T)

E(g)$weight
g<-as.undirected(g)
cl<- cluster_fast_greedy(g ,weights  = E(g)$weight )
#cl<- cluster_walktrap(g ,weights  = E(g)$weight )
#cl<- cluster_louvain(g ,weights  = E(g)$weight )


Membership<-membership(cl)
Membership_fr<-factor(Membership)

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
p2<- p2+ geom_point(aes(col = Membership_fr, shape = Membership_fr), size = 3) 
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


pew_fr<- agg_zeta_fr[agg_zeta_fr$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart"), ]
pew_fr$score_a<- c(3.7, 3.4, 4, 3.3, 2.5, 2.3, 4.1)
pew_fr$score_b<- c(3.2,2.9,3.3,2.9, 2.5,2.4, 3.3)

pew_fr$Membership<- factor(as.numeric(Membership[agg_zeta_fr$names %in% c("TF1", "Le Figaro",  "BFMTV",  "L'Express" , "Le Monde", "Libération",  "Mediapart")]))
pew_fr$Country<- "FR"

#pew_fr$X1<- pew_fr$X1 -  mean(pew_fr$X1)
pew_fr$score_a<- pew_fr$score_a -  mean(pew_fr$score_a)
#########


load("Code/05-Static/RESULT_single_it.RData")
load("Data/Static/Data_Env_single_IT.RData")
pages_names_it<- c(unique(EL_x$i)) # a vector of page names

top_names<-c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" )

el<-EL_princ[,1:3]
colnames(el)<-c("from","to", "weight")
g<-graph_from_data_frame(el)

A<-as.matrix(as_adjacency_matrix(g, attr= "weight"))
A= A+t(A)


E(g)$weight
g<-as.undirected(g)
cl<- cluster_fast_greedy(g ,weights  = E(g)$weight )
#cl<- cluster_walktrap(g ,weights  = E(g)$weight )
#cl<- cluster_louvain(g ,weights  = E(g)$weight )

Membership<-membership(cl)
Membership_it<-factor(Membership)

beta_it<- result[[1]]

zeta_it<- result[[3]]
zeta_it<-data.frame(zeta_it)
colnames(zeta_it) <-c("X1","X2", "i", "it")

zeta_it_sub <-zeta_it[zeta_it$it %in% seq(2000, 15000, 1), ]
zeta_it<-zeta_it[,-2]

agg_zeta_it<-aggregate(zeta_it_sub[,1:2], by = list(zeta_it_sub$i), FUN = mean )
agg_beta_it<-colMeans(beta_it[seq(2000, 15000, 1), ])
agg_zeta_it$names<- pages_names_it

db_it<-data.frame(agg_zeta_it= agg_zeta_it$X1,agg_beta_it= agg_beta_it, names=pages_names_it )
p3<-ggplot(db_it, aes(x = agg_zeta_it, y = agg_beta_it, label = names))
p3<- p3+ geom_point(aes(col = Membership_it, shape = Membership_it), size = 3) 
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


pew_it<- agg_zeta_it[agg_zeta_it$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" ),]
pew_it$score_a<- c(3.3, 3.0, 4.1, 2.7, 3.1, 4.2, 2.9, 4.5)
pew_it$score_b<-c(3.2, 3.2, 3.7, 3.0, 3.2, 3.6, 3.3, 3.7)

#pew_it$X1<- pew_it$X1 -  mean(pew_it$X1)
pew_it$score_a<- pew_it$score_a -  mean(pew_it$score_a)
pew_it$Membership<- factor(as.numeric(Membership[agg_zeta_it$names %in% c("Corriere.della.Sera" , "Il.Fatto.Quotidiano", "Il.Giornale", "la.Repubblica", "La7", "Libero", "Rainews.it", "Tgcom24" )]))
pew_it$Country<- "IT"
##########

load("Code/05-Static/RESULT_single_sp.RData")
load("Data/Static/Data_Env_single_SP.RData")
pages_names_sp<- c(unique(EL_x$i)) # a vector of page names

el<-EL_princ[,1:3]
colnames(el)<-c("from","to", "weight")
g<-graph_from_data_frame(el)

A<-as.matrix(as_adjacency_matrix(g, attr= "weight"))
A= A+t(A)


g<-graph_from_adjacency_matrix(A,weighted = T)

E(g)$weight
g<-as.undirected(g)
cl<- cluster_fast_greedy(g ,weights  = E(g)$weight )
#cl<- cluster_walktrap(g ,weights  = E(g)$weight )
#cl<- cluster_louvain(g ,weights  = E(g)$weight )

Membership<-membership(cl)
Membership_sp<-factor(Membership)

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


top_names<- c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" )


db_sp<-data.frame(agg_zeta_sp= agg_zeta_sp$X1,agg_beta_sp= agg_beta_sp, names=pages_names_sp )
p4<-ggplot(db_sp, aes(x = agg_zeta_sp, y = agg_beta_sp, label = names))
p4<- p4 + geom_point(aes(col = Membership_sp, shape = Membership_sp), size = 3) 
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


#pew_sp$X1<- pew_sp$X1 -  mean(pew_sp$X1)
pew_sp$score_a<- pew_sp$score_a -  mean(pew_sp$score_a)
pew_sp$Membership<- factor(as.numeric(Membership[agg_zeta_sp$names %in% c("ABC.es", "El Mundo", "Antena 3", "La Vanguardia", "El País" )]))
pew_sp$Country<- "SP"
  
PEW_all<-rbind(pew_de, pew_fr, pew_it, pew_sp)

##########

PEW_all$Membership_C<- factor(paste0(as.numeric(PEW_all$Membership),"-",PEW_all$Country))
PEW_all$names_c<- paste0(PEW_all$names, "-", PEW_all$Country)

annotations <- data.frame(
  xpos = c(-Inf,-Inf,Inf,Inf),
  ypos =  c(-Inf, Inf,-Inf,Inf),
  annotateText = c("B-L","T-L"
                   ,"B-R","T-R"),
  hjustvar = c(-.5,   #shifts bottom left 'Text' to the right; make more negative to move it further right
               -.5,   #shifts top left 'tExt' to the right; make more negative to move it further right
               1.5,   #shifts bottom right 'teXt' to the left; make more positive to move it further left
               1.5),  #shifts top right 'texT' to the left; make more positive to move it further left
  vjustvar = c(-1,    #shifts bottom left 'Text' upward; make more negative to move it further up
               2,     #shifts top left 'tExt' downward; make more positive to move it further down
               -1,    #shifts bottom right 'teXt' upward; make more negative to move it further up
               2)     #shifts top right 'texT' downward; make more positive to move it further down
)

p5<- ggplot(PEW_all) +geom_point(aes(x =  X1 , y = score_a, col = Membership_C, shape = Membership_C), size = 4) +geom_text_repel(aes(x =  X1 , y = score_a, label = names_c ), size = 2) + theme_minimal()
p5<- p5 + labs(x = "Latent Leaning", y = "PEW score") # + geom_abline(intercept = 0, slope = 1, linetype = "dashed") 
p5<- p5 + stat_cor(aes(x =  X1 , y = score_a, label = ..r.label..)) #+ xlim(-0.25, 0.25)  
p5<- p5 +scale_shape_manual(values = 0:10)
p5<- p5 + geom_vline(xintercept = 0, linetype = "dashed")
p5<- p5 + geom_hline(yintercept = 0, linetype = "dashed")
p5<- p5 + geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size = 3)
p5 <- p5 + theme_minimal() + theme(legend.position="none", plot.title = element_text(face = "italic",
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


p<- p2 + p1 + p3 + p4 + p5+ plot_annotation(tag_levels = 'A', theme=theme(plot.title=element_text(hjust=0.5))) + plot_layout(design = design)
p

ggsave(p, filename = "Figures/Static/Comparison_FastGreedy.pdf", units = "cm", width = 16*2, height = 9*2 )

