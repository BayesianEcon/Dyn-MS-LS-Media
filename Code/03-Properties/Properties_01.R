
rm(list = ls())

#load the libraries

list.of.packages <- c("ggplot2", "patchwork", "metR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")


library(ggplot2)
library(patchwork)
library(metR)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

#LS property plot

avgDeg = function(n, alpha, d, sigma2_l, sigma2_h, p ){
  
  res_l =(n-1)*exp(alpha)*(4*sigma2_l+1)^(-d/2)
  res_h =(n-1)*exp(alpha)*(4*sigma2_h+1)^(-d/2)
  res = p*res_l + (1-p)*res_h
  
  return(res)
}

n = 100
alpha = seq(0.01, 2.01, 0.01)

sigma2_l = seq(0.01, 3.01, 0.01)
sigma2_h = 4



mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) avgDeg(n, x, d = 1, y,sigma2_h ,1 ) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

breaks = seq(100, 400, by = 100)

require(reshape2)
mmat1 <- melt(mat)
str(mmat1) # to see the names in the melted matrix
g1 <- ggplot(mmat1, aes(x=Var1, y=Var2, z=value) )
g1 <- g1+geom_raster(aes(fill = value)) 
g1 <- g1+stat_contour(aes(col =  after_stat(level)) , col = "black" , lwd = 0.3, breaks = breaks)
g1<-g1+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "Avg. Str.", limits = c(0, 800))
g1 <- g1+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
g1 <- g1+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 1)))
g1 <- g1  +xlim(0,2.05)+ylim(0,3.05) + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                 strip.text.x = element_text(size = 24,face = "bold"),
                                                                                 
                                                                                 axis.title = element_text( size = rel(1)),
                                                                                 axis.title.y = element_text(size = 14),
                                                                                 axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                 axis.text=element_text(size=12),
                                                                                 axis.line = element_blank(),
                                                                                 axis.ticks = element_line(),
                                                                                 panel.grid = element_blank(),
                                                                                 legend.title = element_text(face="italic", size=12),
                                                                                 legend.text = element_text(size=8),
                                                                                 panel.border = element_rect(colour = "black", fill = NA))








g1

n = 100
alpha = seq(0.01, 2.01, 0.01)

sigma2_l = seq(0.01, 3.01, 0.01)
sigma2_h = 4



mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) avgDeg(n, x, d = 1, y,sigma2_h ,0.75 ) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

require(reshape2)
mmat1 <- melt(mat)
str(mmat1) # to see the names in the melted matrix
g2 <- ggplot(mmat1, aes(x=Var1, y=Var2, z=value) )
g2 <- g2+geom_raster(aes(fill = value)) 
g2 <- g2+stat_contour(aes(col =  after_stat(level)) , col = "black" , lwd = 0.3,  breaks = breaks)
g2<-g2+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "Avg. Str.", limits = c(0, 800))
g2 <- g2+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
g2 <- g2+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 0.75)))
g2 <- g2  +xlim(0,2.05)+ylim(0,3.05) + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                 strip.text.x = element_text(size = 24,face = "bold"),
                                                                                 
                                                                                 axis.title = element_text( size = rel(1)),
                                                                                 axis.title.y = element_text(size = 14),
                                                                                 axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                 axis.text=element_text(size=12),
                                                                                 axis.line = element_blank(),
                                                                                 axis.ticks = element_line(),
                                                                                 panel.grid = element_blank(),
                                                                                 legend.title = element_text(face="italic", size=12),
                                                                                 legend.text = element_text(size=8),
                                                                                 panel.border = element_rect(colour = "black", fill = NA))







g2



n = 100
alpha = seq(0.01, 2.01, 0.01)

sigma2_l = seq(0.01, 3.01, 0.01)
sigma2_h = 4



mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) avgDeg(n, x, d = 1, y,sigma2_h ,0.5 ) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

require(reshape2)
mmat1 <- melt(mat)
str(mmat1) # to see the names in the melted matrix
g3 <- ggplot(mmat1, aes(x=Var1, y=Var2, z=value) )
g3 <- g3+geom_raster(aes(fill = value)) 
g3 <- g3+stat_contour(aes(col =  after_stat(level)) , col = "black" , lwd = 0.3, breaks = breaks)
g3<-g3+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "Avg. Str.", limits = c(0, 800))
g3 <- g3+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
g3 <- g3+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste(q[L]== 0.5)))
g3 <- g3  +xlim(0,2.05)+ylim(0,3.05) + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                               strip.text.x = element_text(size = 24,face = "bold"),
                                                               
                                                               axis.title = element_text( size = rel(1)),
                                                               axis.title.y = element_text(size = 14),
                                                               axis.title.x = element_text(size = 14, vjust = -0.2),
                                                               axis.text=element_text(size=12),
                                                               axis.line = element_blank(),
                                                               axis.ticks = element_line(),
                                                               panel.grid = element_blank(),
                                                               legend.title = element_text(face="italic", size=12),
                                                               legend.text = element_text(size=8),
                                                               panel.border = element_rect(colour = "black", fill = NA))





g3


#### Dispersion Index 

alpha = seq(0.01, 2.01, 0.01)
sigma2 = seq(0.01, 3.01, 0.01)



DI = function(n, alpha, d, sigma2_l, sigma2_h, p){
  
  #a = 2*(n-1)*exp(2*alpha)*(2*sigma2 + 0.5 )^(-d/2)*( (4/(4+sigma2)) + (1/sigma2) )^(-d/2)*(sigma2)^(-d/2)*2^(-d/2) +  (n-1)*(n-2)*exp(2*alpha)*(sigma2+0.5)^(-d)*( (2/(2+sigma2)) + (1/sigma2))^(-d/2)*(sigma2)^(-d/2)*2^(-d)
  
  a_l = exp(2*alpha)*(n-1)*(8*sigma2_l+1)^(-d/2)+(n-1)*(n-2)*exp(2*alpha)*(2*sigma2_l+1)^(-d/2)*(6*sigma2_l+1)^(-d/2)
  a_h = exp(2*alpha)*(n-1)*(8*sigma2_h+1)^(-d/2)+(n-1)*(n-2)*exp(2*alpha)*(2*sigma2_h+1)^(-d/2)*(6*sigma2_h+1)^(-d/2)
  
  a = p*a_l+(1-p)*a_h
  
  b_l = (n-1)*exp(alpha)*(4*sigma2_l+1)^(-d/2)
  b_h = (n-1)*exp(alpha)*(4*sigma2_h+1)^(-d/2)
  
  b = p*b_l+(1-p)*b_h
  
  di =  1+ (a/b) - b
  
  return(di)
  
}




mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) DI(n=100, x, d = 1, y, sigma2_h ,1) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

breaks = seq(10, 25, by = 5)

require(reshape2)
mmat4 <- melt(mat)
q1 <- ggplot(mmat4, aes(x=Var1, y=Var2, z=value) )
q1 <- q1 +geom_raster(aes(fill = value)) 
q1 <- q1+stat_contour(aes(col = after_stat(level)), col = "black" , lwd = 0.3, breaks = breaks)
q1<-  q1+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "D.I.", limits = c(0, 300))
q1 <- q1+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 90))
q1 <- q1+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 1)))
q1 <- q1 + xlim(0,2.05)+ylim(0,3.05) + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                               strip.text.x = element_text(size = 24,face = "bold"),
                                                               
                                                               axis.title = element_text( size = rel(1)),
                                                               axis.title.y = element_text(size = 14),
                                                               axis.title.x = element_text(size = 14, vjust = -0.2),
                                                               axis.text=element_text(size=12),
                                                               axis.line = element_blank(),
                                                               axis.ticks = element_line(),
                                                               panel.grid = element_blank(),
                                                               legend.title = element_text(face="italic", size=12),
                                                               legend.text = element_text(size=8),
                                                               panel.border = element_rect(colour = "black", fill = NA))







q1


mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) DI(n=100, x, d = 1, y, sigma2_h , 0.75) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

breaks = seq(20, 80, by = 20)


require(reshape2)
mmat4 <- melt(mat)
q2 <- ggplot(mmat4, aes(x=Var1, y=Var2, z=value) )
q2 <- q2 +geom_raster(aes(fill = value)) 
q2 <- q2+stat_contour(aes(col = ..level..), col = "black" , lwd = 0.3,  breaks = breaks)
q2<-  q2+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "D.I.", limits = c(0, 300))
q2 <- q2+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
q2 <- q2+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 0.75)))
q2 <- q2 + xlim(0,2.05)+ylim(0,3.05) + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                 strip.text.x = element_text(size = 24,face = "bold"),
                                                                                 
                                                                                 axis.title = element_text( size = rel(1)),
                                                                                 axis.title.y = element_text(size = 14),
                                                                                 axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                 axis.text=element_text(size=12),
                                                                                 axis.line = element_blank(),
                                                                                 axis.ticks = element_line(),
                                                                                 panel.grid = element_blank(),
                                                                                 legend.title = element_text(face="italic", size=12),
                                                                                 legend.text = element_text(size=8),
                                                                                 panel.border = element_rect(colour = "black", fill = NA))






q2



mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) DI(n=100, x, d = 1, y, sigma2_h ,0.25) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

breaks = seq(40, 100, by = 20)


require(reshape2)
mmat4 <- melt(mat)
q3 <- ggplot(mmat4, aes(x=Var1, y=Var2, z=value) )
q3 <- q3 +geom_raster(aes(fill = value)) 
q3 <- q3+stat_contour(aes(col = ..level..), col = "black" , lwd = 0.3,breaks=breaks)
q3<-  q3+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "D.I.", limits = c(0, 300))
q3 <- q3+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
q3 <- q3+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 0.5)))
q3 <- q3 + xlim(0,2.05)+ylim(0,3.05) + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                 strip.text.x = element_text(size = 24,face = "bold"),
                                                                                 
                                                                                 axis.title = element_text( size = rel(1)),
                                                                                 axis.title.y = element_text(size = 14),
                                                                                 axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                 axis.text=element_text(size=12),
                                                                                 axis.line = element_blank(),
                                                                                 axis.ticks = element_line(),
                                                                                 panel.grid = element_blank(),
                                                                                 legend.title = element_text(face="italic", size=12),
                                                                                 legend.text = element_text(size=8),
                                                                                 panel.border = element_rect(colour = "black", fill = NA))





q3


Var = function(n, alpha, d, sigma2_l, sigma2_h , p){
  
  a_l = exp(2*alpha)*(n-1)*(8*sigma2_l+1)^(-d/2)+(n-1)*(n-2)*exp(2*alpha)*(2*sigma2_l+1)^(-d/2)*(6*sigma2_l+1)^(-d/2)
  a_h = exp(2*alpha)*(n-1)*(8*sigma2_h+1)^(-d/2)+(n-1)*(n-2)*exp(2*alpha)*(2*sigma2_h+1)^(-d/2)*(6*sigma2_h+1)^(-d/2)
  
  a = p*a_l+(1-p)*a_h
  
  b_l = (n-1)*exp(alpha)*(4*sigma2_l+1)^(-d/2)
  b_h = (n-1)*exp(alpha)*(4*sigma2_h+1)^(-d/2)
  
  b = p*b_l+(1-p)*b_h
  
  di =  a +b - b^2
  
  return(sqrt(di))
  
}



mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) Var(n=100, x, d = 1, y, sigma2_h, 1) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

breaks = seq(40, 100, 20)

require(reshape2)
mmat3 <- melt(mat)
p1 <- ggplot(mmat3, aes(x=Var1, y=Var2, z=value) )
p1 <- p1+geom_raster(aes(fill = value)) 
p1 <- p1+stat_contour(aes(col = ..level..), col = "black" , lwd = 0.3, breaks = breaks)
p1<-  p1+ scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "St.D.", limits = c(0, 400))
p1 <- p1+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
p1 <- p1  + labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 1)))
p1 <- p1 + xlim(0,2.05)+ylim(0,3.05)  + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                  strip.text.x = element_text(size = 24,face = "bold"),
                                                                                  
                                                                                  axis.title = element_text( size = rel(1)),
                                                                                  axis.title.y = element_text(size = 14),
                                                                                  axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                  axis.text=element_text(size=12),
                                                                                  axis.line = element_blank(),
                                                                                  axis.ticks = element_line(),
                                                                                  panel.grid = element_blank(),
                                                                                  legend.title = element_text(face="italic", size=12),
                                                                                  legend.text = element_text(size=8),
                                                                                  panel.border = element_rect(colour = "black", fill = NA))





p1

mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) Var(n=100, x, d = 1, y, sigma2_h, 0.75) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

require(lattice)

breaks = seq(60, 120, 20)


require(reshape2)
mmat3 <- melt(mat)
p2  <- ggplot(mmat3, aes(x=Var1, y=Var2, z=value) )
p2  <- p2 +geom_raster(aes(fill = value)) 
p2  <- p2 +stat_contour(aes(col = ..level..), col = "black" , lwd = 0.3, breaks = breaks)
p2 <-  p2 + scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "St.D.", limits = c(0, 400))
p2 <- p2+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
p2 <- p2+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 0.75)))
p2 <- p2 + xlim(0,2.05)+ylim(0,3.05)  + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                  strip.text.x = element_text(size = 24,face = "bold"),
                                                                                  
                                                                                  axis.title = element_text( size = rel(1)),
                                                                                  axis.title.y = element_text(size = 14),
                                                                                  axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                  axis.text=element_text(size=12),
                                                                                  axis.line = element_blank(),
                                                                                  axis.ticks = element_line(),
                                                                                  panel.grid = element_blank(),
                                                                                  legend.title = element_text(face="italic", size=12),
                                                                                  legend.text = element_text(size=8),
                                                                                  panel.border = element_rect(colour = "black", fill = NA))





p2


mat <- outer( alpha,  
              sigma2_l, 
              Vectorize( function(x,y) Var(n=100, x, d = 1, y, sigma2_h, 0.5) ) )

rownames(mat)<-alpha
colnames(mat)<-sigma2_l

breaks = seq(80, 140, 20)



require(reshape2)
mmat3 <- melt(mat)
p3  <- ggplot(mmat3, aes(x=Var1, y=Var2, z=value) )
p3  <- p3 +geom_raster(aes(fill = value)) 
p3  <- p3 +stat_contour(aes(col = ..level..), col = "black" , lwd = 0.3,  breaks = breaks)
p3 <-  p3 + scale_fill_gradientn(colours=c("#fcfbfa","#da6d50"), name = "St.D.", limits = c(0, 400))
p3 <- p3+geom_text_contour(check_overlap = TRUE, breaks = breaks,  label.placer = label_placer_flattest(ref_angle = 45))
p3 <- p3+ labs(x = expression(alpha), y = expression(sigma[L]^{2}*beta), title = expression(paste( q[L]== 0.5)))
p3 <- p3 + xlim(0,2.05)+ylim(0,3.05)  + theme_minimal() + theme_minimal() + theme(legend.position="bottom", plot.title = element_text( size = 14),
                                                                                  strip.text.x = element_text(size = 24,face = "bold"),
                                                                                  
                                                                                  axis.title = element_text( size = rel(1)),
                                                                                  axis.title.y = element_text(size = 14),
                                                                                  axis.title.x = element_text(size = 14, vjust = -0.2),
                                                                                  axis.text=element_text(size=12),
                                                                                  axis.line = element_blank(),
                                                                                  axis.ticks = element_line(),
                                                                                  panel.grid = element_blank(),
                                                                                  legend.title = element_text(face="italic", size=12),
                                                                                  legend.text = element_text(size=8),
                                                                                  panel.border = element_rect(colour = "black", fill = NA))






p3



g1 <- (g1|g2|g3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
g2<-(p1|p2|p3) + plot_layout(guides = "collect")   & theme(legend.position = "bottom")
g3<-(q1|q2|q3) + plot_layout(guides = "collect")   & theme(legend.position = "bottom")
g<-(g1/g2/g3)+ plot_annotation(tag_levels = "A") & theme(plot.title = element_text( size = 12)) #+ plot_layout(guides = "collect") &  scale_fill_gradientn(limits = range(c(mmat1$value, mmat3$value, mmat4$value)) , colors=c("yellow","red"), name = "Avg. Degree") & theme(legend.position = "bottom")

g

ggsave(g, file = "Figures/Properties/Figure4.pdf", width = 21, height = 29.7, unit = "cm")



