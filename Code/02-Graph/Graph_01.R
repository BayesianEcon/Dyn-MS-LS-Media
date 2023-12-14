 
rm(list = ls())

#load the libraries

list.of.packages <- c("lubridate", "plyr", "data.table", "ggplot2", "igraph", "dplyr", "scales", "grid",
                      "ggraph","ggiraphExtra","ggtext","ggnewscale","tidygraph", "ggpmisc", "patchwork")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")



library(lubridate)
library(plyr)
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(grid)
library(igraph)
library(ggraph)
library(ggiraphExtra)
library(ggtext)
library(ggnewscale)
library(tidygraph)
library(ggpmisc)
library(patchwork)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################


load("Data/Dynamic/DataEnv_IT_all.RData")

#################################################

DBplane$date<-as_date(DBplane$t)
DBplane$date<-ymd(DBplane$date)
DBplane$y<-year(DBplane$date)
DBplane$m<-month(DBplane$date)
DBplane$m<-as.character(DBplane$m)

sub0 = function(x){
  
  if(nchar(x)== 1){
    
    x = paste0("0",x)
    
  }else{x}
}

DBplane$m<- unlist(sapply(DBplane$m, FUN = sub0))

DBplane$ym<- paste0(DBplane$y,"_",DBplane$m)
DBplane$ym<- as.factor(DBplane$ym)
DBplane$leaning<- DBplane$leaning - mean(DBplane$leaning)


date_one<- c("2016-7", "2016-8")
date_two<- c("2016-9", "2016-10")

date_one<- "2015"
date_two<- "2016"

DBplane_2015 <- DBplane[DBplane$y %in% date_one,]
DBplane_2016 <- DBplane[DBplane$y %in% date_two,]

EL_princ$date<-as_date(EL_princ$tindex)
EL_princ$y<-year(EL_princ$date)
EL_princ$m<-month(EL_princ$date)

EL_princ_2015<-EL_princ[EL_princ$y %in% date_one,]
EL_princ_2016<-EL_princ[EL_princ$y %in% date_two,]


#########################

EL_princ_2015_agg<- aggregate(EL_princ_2015$w, by = list(i = EL_princ_2015$i, EL_princ_2015$j), FUN = sum)

g2015<-graph_from_edgelist(as.matrix(EL_princ_2015_agg[, 1:2]), directed = FALSE)
E(g2015)$weight<-EL_princ_2015_agg$x +1
leaning<- aggregate(DBplane_2015$leaning, by = list(name =DBplane_2015$fb_name ), FUN = mean)[,2]

leaning2015 <-  leaning
V(g2015)$leaning<- leaning2015
V(g2015)$leaning<-(V(g2015)$leaning - min(V(g2015)$leaning ))/(max(V(g2015)$leaning) - min(V(g2015)$leaning))
leaning2015<-V(g2015)$leaning


V(g2015)$cl = factor(membership(cluster_fast_greedy(g2015)))
V(g2015)$ncomments<- aggregate(DBplane_2015$ncomments, by = list(name =DBplane_2015$fb_name ), FUN = sum)[,2]

adj.matrix = get.adjacency(g2015,attr = "weight", sparse=FALSE)
adj.matrix = t(apply(adj.matrix,1, FUN = function(x){x/sum(x)}))

km <- kmeans(adj.matrix , centers = 2)


V(g2015)$sign<- leaning2015

score<-t(adj.matrix)%*%V(g2015)$sign
score<- (score - min(score))/(max(score) - min(score))
score<- score - mean(score)


score_save<-rep(0, 1000)


# if(score[370){score<- -1*score}

score_save<-rep(0, 1000)

q<- 0.8
for(i in 1:1000){
  score<- q*t(adj.matrix)%*%score + (1-q)*leaning2015
  score<- (score - min(score))/(max(score) - min(score))
  score<- score - mean(score)
  
  
  score_save[i]<-score[10]
  print(i)
}

score2015<-  score
V(g2015)$score <- score2015

cl<- fastgreedy.community(g2015 ,weights  = E(g2015)$weight )



weights <- ifelse(igraph::crossing(cl, g2015), 1, 4)
ly <- layout_with_fr(g2015, weights=weights)

or_coords= data.frame(ly[,1:2])

V(g2015)$cx<-or_coords[,1]
V(g2015)$cy<-or_coords[,2]

V(g2015)$name<-gsub("\\."," ", V(g2015)$name)


############################

EL_princ_2016_agg<- aggregate(EL_princ_2016$w, by = list(i = EL_princ_2016$i, EL_princ_2016$j), FUN = sum)

g2016<-graph_from_edgelist(as.matrix(EL_princ_2016_agg[, 1:2]), directed = FALSE)
E(g2016)$weight<-EL_princ_2016_agg$x +1
leaning<- aggregate(DBplane_2016$leaning, by = list(name = DBplane_2016$fb_name ), FUN = mean)[,2]
leaning2016 <- leaning

V(g2016)$leaning<-leaning2016
V(g2016)$leaning<- (V(g2016)$leaning - min(V(g2016)$leaning ))/(max(V(g2016)$leaning) - min(V(g2016)$leaning))
leaning2016<-V(g2016)$leaning

V(g2016)$cl = factor(membership(cluster_fast_greedy(g2016)))
V(g2016)$ncomments<- aggregate(DBplane_2016$ncomments, by = list(name =DBplane_2016$fb_name ), FUN = sum)[,2]


adj.matrix = get.adjacency(g2016,attr = "weight", sparse=FALSE)
adj.matrix = t(apply(adj.matrix,1, FUN = function(x){x/sum(x)}))


km <- kmeans(adj.matrix , centers = 2)


V(g2016)$sign<- leaning2016
V(g2016)$name<-gsub("\\."," ", V(g2016)$name)


score<-t(adj.matrix)%*%V(g2016)$sign
score<- (score - min(score))/(max(score) - min(score))
score<- score - mean(score)


#if(score[37]>0){score<- -1*score}

score_save<-rep(0, 1000)

q<- 0.8
for(i in 1:1000){
  score<- q*t(adj.matrix)%*%score + (1-q)*leaning2016
  score<- (score - min(score))/(max(score) - min(score))
  score<- score - mean(score)
  
  
  score_save[i]<-score[10]
  print(i)
}

score2016<-  score

# p1 <- hist(leaning2015, breaks= 10)                     # centered at 4
# p2 <- hist(leaning2016, breaks= 10)                     # centered at 6
# plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1))  # first histogram
# plot( p2, col=rgb(1,0,0,1/4), add=T)  # second

# cl<- make_clusters(
#   g2016,
#   membership = km$cluster
# )
# 
# weights <- 1/c(dist(leaning2016))^2
# 
# ly<-create_layout(g2016, "fr",  weights=weights)

V(g2016)$score<-score2016
V(g2016)$cx<-or_coords$X1
V(g2016)$cy<-or_coords$X2

####################################

V(g2015)$norm_ncomments<- (V(g2015)$ncomments - min( c(V(g2015)$ncomments, V(g2016)$ncomments) ))/(max(c(V(g2015)$ncomments, V(g2016)$ncomments)) - min( c(V(g2015)$ncomments, V(g2016)$ncomments) ))
V(g2016)$norm_ncomments<- (V(g2016)$ncomments - min( c(V(g2015)$ncomments, V(g2016)$ncomments) ))/(max(c(V(g2015)$ncomments, V(g2016)$ncomments)) - min( c(V(g2015)$ncomments, V(g2016)$ncomments) ))

V(g2015)$norm_score<- (V(g2015)$score - min( c(V(g2015)$score, V(g2016)$score) ))/(max(c(V(g2015)$score, V(g2016)$score)) - min( c(V(g2015)$score, V(g2016)$score) ))
V(g2016)$norm_score<- (V(g2016)$score - min( c(V(g2015)$score, V(g2016)$score) ))/(max(c(V(g2015)$score, V(g2016)$score)) - min( c(V(g2015)$score, V(g2016)$score) ))

E(g2015)$norm_weight<- (E(g2015)$weight - min( c(E(g2015)$weight, E(g2016)$weight) ))/(max(c(E(g2015)$weight, E(g2016)$weight)) - min( c(E(g2015)$weight, E(g2016)$weight) ))
E(g2016)$norm_weight<- (E(g2016)$weight - min( c(E(g2015)$weight, E(g2016)$weight) ))/(max(c(E(g2015)$weight, E(g2016)$weight)) - min( c(E(g2015)$weight, E(g2016)$weight) ))

#  leaning2015<- (leaning2015 - min( c(leaning2015, leaning2016) ))/(max(c(leaning2015,leaning2016)) - min( c(leaning2015, leaning2016) ))
# leaning2016<- (leaning2016 - min( c(leaning2015, leaning2016) ))/(max(c(leaning2015,leaning2016)) - min( c(leaning2015, leaning2016) ))
leaning2015<- (leaning2015 - min( c(leaning2015) ))/(max(c(leaning2015)) - min( c(leaning2015) ))
leaning2016<- (leaning2016 - min( c(leaning2016) ))/(max(c(leaning2016)) - min( c(leaning2016) ))


palette=c("#D9D9DA","#DC143C","#1446A0")

####################################
coords_2015<-data.frame(x = V(g2015)$cx, y = V(g2015)$cy)
tbl2015 = as_tbl_graph(g2015)



gr_all=ggraph(tbl2015,  coords_2015) + 
  geom_edge_arc(aes(alpha=norm_weight),strength = 0.15, show.legend = F) + 
  scale_edge_alpha_continuous(range = c(0.001,2) , limits = c(0,1))+
  geom_node_point(aes(fill=  leaning2015, size= norm_ncomments) , shape = 21)+
  #scale_fill_manual(values = c("red", "blue"))+
  #geom_hline(yintercept = 1.4, col = "white")+
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",   midpoint = 0.5)+
  #scale_fill_gradient2(low = "red", mid = "white", high = "blue",   midpoint = 0.5, limits = c(0,1))+
  geom_node_text(aes(label = name, size = norm_ncomments), repel = T, fontface = "italic", show.legend = F)+
  scale_size_continuous(limits = c(0,1))+
  labs(title = paste("Year: " , date_one), subtitle = "Audience-Duplication Network",
       size= expression(italic("Engagement")), fill=expression(italic("Network-Based Leaning")))+
  theme_void()+
  theme(legend.position = "bottom",
        text=element_text(size=18),
        title = element_text(hjust = 0.5, face = "italic"),
        plot.subtitle  = element_text( face = "italic"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.box="vertical", legend.margin=margin()
  )+
  #coord_equal()+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5),
         size=guide_legend(title.position = "top", title.hjust = 0.5))
#size = guide_legend(override.aes = list(size = 5)))
#geom_node_point(aes(color=factor(reliability)))
gr_all


pp2015<-gr_all



#####################


coords_2016<-data.frame(x = V(g2016)$cx, y = V(g2016)$cy )

tbl2016 = as_tbl_graph(g2016)
#devtools::install_github('thomasp85/ggraph')

gr_all=ggraph(tbl2016,  coords_2016) + 
  geom_edge_arc(aes(alpha=norm_weight),strength = 0.15, show.legend = F) + 
  scale_edge_alpha_continuous(range = c(0.001,2), limits = c(0,1))+
  geom_node_point(aes(fill=  leaning2016, size= norm_ncomments) , shape = 21)+
  #geom_hline(yintercept = 1.4, col = "white")+
  #scale_fill_manual(values = c("blue", "red"))+
  scale_fill_gradient2(low = "red", mid = "white", high = "blue",   midpoint = 0.5)+
  # scale_fill_gradient2(low = "red", mid = "white", high = "blue",   midpoint = 0.5, limits= c(0,1))+
  geom_node_text(aes(label = name, size = norm_ncomments), repel = T, fontface = "italic", show.legend = F)+
  scale_size_continuous(limits = c(0,1))+
  labs(title =   paste("Year: " , date_two), subtitle = "Audience-Duplication Network",
       size= expression(italic("Engagement")), fill=expression(italic("Network-Based Leaning")))+
  theme_void()+
  theme(legend.position = "bottom",
        text=element_text(size=18),
        title = element_text(hjust = 0.5, face = "italic"),
        plot.subtitle  = element_text( face = "italic"),
        legend.title = element_text(size=18),
        legend.text = element_text(size=14),
        legend.box="vertical", legend.margin=margin()
  )+
  #coord_equal()+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5),
         size=guide_legend(title.position = "top", title.hjust = 0.5))
#size = guide_legend(override.aes = list(size = 5)))
#geom_node_point(aes(color=factor(reliability)))
gr_all


pp2016<-gr_all


###############

EL_princ$date<-as_date(EL_princ$tindex)
EL_princ$y<-year(EL_princ$date)
EL_princ$m<-month(EL_princ$date)
EL_princ$m<- unlist(sapply(EL_princ$m, FUN = sub0))
EL_princ$ym<- paste0(EL_princ$y,"-",EL_princ$m)

count_list<-aggregate(EL_princ$w, by = list(date = EL_princ$date), FUN = mean)
count_tot<-aggregate(DBplane$ncomments, by = list(date = DBplane$date), FUN = sum)

#count_list$x<-  count_list$x/count_tot$x

count_list$year<-year(count_list$date)
count_list$m<-month(count_list$date)
count_list$m<- unlist(sapply(count_list$m, FUN = sub0))
count_list$ym<- factor(paste0(count_list$y,"-",count_list$m))

tabs<-data.frame(aggregate(count_list$x, by =  list(ym=count_list$ym), FUN = mean))

my_custom_labels <- as.character(tabs$ym)
my_custom_labels[seq(2, 24,2)]<- ""

box<-ggplot(count_list, aes(x = ym, y = x))+stat_summary( width=0.4, col = "red",fun.min = function(z) { quantile(z,0.25) },
                                                          fun.max = function(z) { quantile(z,0.75) },
                                                          fun = median, geom = "errorbar") 
box<-box + geom_point(data=tabs, aes(x = factor(ym), y =x, group = 1 ), color = "black", linetype = "dashed") 
box<- box + labs(title = "", subtitle = "Distribution of Average Weight over Time", x = "Text-based Leaning")
box<- box + scale_x_discrete(labels = my_custom_labels)  + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( face = "italic", size = 16),
                                                                                   plot.subtitle = element_text( face = "italic", size = 12),
                                                                                   
                                                                                   strip.text.x = element_text(size = 12,face = "bold"),
                                                                                   
                                                                                   axis.title = element_text( face = "italic",size = 12),
                                                                                   axis.title.y = element_blank(),
                                                                                   #axis.title.x = element_text(size = 18, vjust = -0.2),
                                                                                   axis.title.x=element_blank(),
                                                                                   axis.text.x=element_text( size = 16 ),
                                                                                   axis.text.y=element_text( size = 16),
                                                                                   axis.line = element_blank(),
                                                                                   axis.ticks = element_line(),
                                                                                   panel.grid = element_blank(),
                                                                                   legend.title = element_text(face="italic", size=12),
                                                                                   legend.text = element_text(size=8))

box

layout2 <- "AAAABBBB
            AAAABBBB
            AAAABBBB
            CCCCCCCC
            CCCCCCCC"


resa<- pp2015 + pp2016 + box + plot_layout(guides = "collect", design= layout2) & theme(legend.position = "none", axis.title.y = element_blank(), plot.title = element_text(size = 20), plot.subtitle = element_text(size = 16)  )
resa


ggsave(resa, file ="Figures/Graph/Figure1.pdf", unit = "cm", width = 16*2, height = 9*2)
