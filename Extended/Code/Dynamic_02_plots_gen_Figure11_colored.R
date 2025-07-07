library(mvtnorm)
library(Rcpp)
library(RcppDist)
library(RcppParallel)
library(RcppArmadillo)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(scales)
library(dplyr)
library(data.table)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(lubridate)

load("/Extended/Data/dates.RData")
load("/Extended/Data/plot_df_fr.RData")
load("/Extended/Data/xi_fr_agg_fr.RData")
load("/Extended/Data/EY_fr.RData")

 K= 5

q <- ggplot(plot_df, aes(x = Time, y = Value, group = Time, color = State)) +
  geom_boxplot(
    size = 0.2, fill = NA, outlier.size = 0.3,
    outlier.shape = 16, outlier.alpha = 0.1
  ) +
  stat_summary(
    fun = mean, geom = "step", aes(group = 1),
    size = 0.1, color = "gray70", linetype = "solid"
  ) +
  geom_point(
    data = data.frame(Time = dates, Value = EY),
    aes(x = Time, y = Value),
    inherit.aes = FALSE, size = 0.2, color = "gray30" #, linetype = "solid"
  ) +
  labs(y = "Expected Strength", title = "France") +
  theme_classic() +
  scale_color_gradientn(colors = c("#0000ff", "#ff0080", "#c20824"  )) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "italic", size = 24, hjust = 0.5),
    strip.text.x = element_text(size = 16, face = "italic"),
    axis.title = element_text(face = "italic", size = 16),
    axis.title.y = element_text(size = 16, vjust = -0.2),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0, 0, -2, 0), "cm"),
    axis.line = element_blank(),
    axis.ticks.length.x = unit(0, 'lines'),
    panel.grid = element_blank(),
    legend.title = element_text(face = "italic", size = 16),
    legend.text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA)
  )
q1<-q

q<-ggplot(xi_fr_agg, aes(x = Time, y = (K+1)-State), col ="red") 
q <- q +  geom_step(size = 0.2, color = "grey70", alpha = 0.8)+ theme_classic()
q <- q +  geom_point(aes(color =  (K+1)-State), size = 0.1, alpha = 0.8)
q<- q+ labs(y = "States") + theme_classic()
q <- q+ scale_color_gradientn(colors=c( "#0000ff", "#ff0080", "#c20824"))
q<- q + theme(legend.position="none", 
              plot.title = element_blank(),
              strip.text.x = element_text(size = 16,face = "italic"),
              axis.title = element_blank(),
              axis.title.y = element_text(size = 16, vjust = -0.2),
              axis.title.x = element_blank(),
              axis.text=element_text(size=5),
              axis.line = element_blank(),
              axis.ticks = element_line(),
              panel.grid = element_blank(),
              #plot.margin = unit(c(-20,0,0,0), "cm"),
              legend.title = element_text(face="italic", size=16),
              legend.text = element_text(size=14),
              panel.border = element_rect(colour = "black", fill = NA))
q2<-q

gv_fr <- (q1 / q2)  +  plot_layout(heights = c(4, 1)) #& theme(plot.margin = unit(c(0,0,-0.1,0), "cm"),)

save(EY , file = "EY_fr.RData")
##########################

load("/Extended/Data/plot_df_de.RData")
load("/Extended/Data/xi_fr_agg_de.RData")
load("/Extended/Data/EY_de.RData")

q <- ggplot(plot_df, aes(x = Time, y = Value, group = Time, color = State)) +
  geom_boxplot(
    size = 0.2, fill = NA, outlier.size = 0.3,
    outlier.shape = 16, outlier.alpha = 0.1
  ) +
  stat_summary(
    fun = mean, geom = "step", aes(group = 1),
    size = 0.1, color = "gray70", linetype = "solid"
  ) +
  geom_point(
    data = data.frame(Time = dates, Value = EY),
    aes(x = Time, y = Value),
    inherit.aes = FALSE, size = 0.2, color = "gray30"#, linetype = "solid"
  ) +
  labs(y = "Expected Strength", title = "Germany") +
  theme_classic() +
  scale_color_gradientn(colors = c("#0000ff",  "#ff0080", "#c20824"  )) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "italic", size = 24, hjust = 0.5),
    strip.text.x = element_text(size = 16, face = "italic"),
    axis.title = element_text(face = "italic", size = 16),
    axis.title.y = element_text(size = 16, vjust = -0.2),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0, 0, -2, 0), "cm"),
    axis.line = element_blank(),
    axis.ticks.length.x = unit(0, 'lines'),
    panel.grid = element_blank(),
    legend.title = element_text(face = "italic", size = 16),
    legend.text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA)
  )
q1<-q

q<-ggplot(xi_fr_agg, aes(x = Time, y = (K+1)-State), col ="red") 
q <- q +  geom_step(size = 0.2, color = "grey70", alpha = 0.8)+ theme_classic()
q <- q +  geom_point(aes(color =  (K+1)-State), size = 0.1, alpha = 0.8)
q<- q+ labs(y = "States") + theme_classic()
q <- q+ scale_color_gradientn(colors=c( "#0000ff", "#ff0080", "#c20824"))
q<- q + theme(legend.position="none", 
              plot.title = element_blank(),
              strip.text.x = element_text(size = 16,face = "italic"),
              axis.title = element_blank(),
              axis.title.y = element_text(size = 16, vjust = -0.2),
              axis.title.x = element_blank(),
              axis.text=element_text(size=5),
              axis.line = element_blank(),
              axis.ticks = element_line(),
              panel.grid = element_blank(),
              #plot.margin = unit(c(-20,0,0,0), "cm"),
              legend.title = element_text(face="italic", size=16),
              legend.text = element_text(size=14),
              panel.border = element_rect(colour = "black", fill = NA))
q2<-q
gv_de <- (q1 / q2)  +  plot_layout(heights = c(4, 1)) #& theme(plot.margin = unit(c(0,0,-0.1,0), "cm"),)



#####################################
load("/Extended/Data/plot_df_it.RData")
load("/Extended/Data/xi_fr_agg_it.RData")
load("/Extended/Data/EY_it.RData")

q <- ggplot(plot_df, aes(x = Time, y = Value, group = Time, color = State)) +
  geom_boxplot(
    size = 0.2, fill = NA, outlier.size = 0.3,
    outlier.shape = 16, outlier.alpha = 0.1
  ) +
  stat_summary(
    fun = mean, geom = "step", aes(group = 1),
    size = 0.1, color = "gray70", linetype = "solid"
  ) +
  geom_point(
    data = data.frame(Time = dates, Value = EY),
    aes(x = Time, y = Value),
    inherit.aes = FALSE, size = 0.2, color = "gray30"#, linetype = "solid"
  ) +
  labs(y = "Expected Strength", title = "Italy") +
  theme_classic() +
  scale_color_gradientn(colors = c("#c20824" , "#ff0080", "#0000ff")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "italic", size = 24, hjust = 0.5),
    strip.text.x = element_text(size = 16, face = "italic"),
    axis.title = element_text(face = "italic", size = 16),
    axis.title.y = element_text(size = 16, vjust = -0.2),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0, 0, -2, 0), "cm"),
    axis.line = element_blank(),
    axis.ticks.length.x = unit(0, 'lines'),
    panel.grid = element_blank(),
    legend.title = element_text(face = "italic", size = 16),
    legend.text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA)
  )
q1<-q

q<-ggplot(xi_fr_agg, aes(x = Time, y = (K+1)-State), col ="red") 
q <- q +  geom_step(size = 0.2, color = "grey70", alpha = 0.8)+ theme_classic()
q <- q +  geom_point(aes(color =  (K+1)-State), size = 0.1, alpha = 0.8)
q<- q+ labs(y = "States") + theme_classic()
q <- q+ scale_color_gradientn(colors=c( "#0000ff", "#ff0080", "#c20824"))
q<- q + theme(legend.position="none", 
              plot.title = element_blank(),
              strip.text.x = element_text(size = 16,face = "italic"),
              axis.title = element_blank(),
              axis.title.y = element_text(size = 16, vjust = -0.2),
              axis.title.x = element_blank(),
              axis.text=element_text(size=5),
              axis.line = element_blank(),
              axis.ticks = element_line(),
              panel.grid = element_blank(),
              #plot.margin = unit(c(-20,0,0,0), "cm"),
              legend.title = element_text(face="italic", size=16),
              legend.text = element_text(size=14),
              panel.border = element_rect(colour = "black", fill = NA))
q2<-q
gv_it <- (q1 / q2)  +  plot_layout(heights = c(4, 1)) #& theme(plot.margin = unit(c(0,0,-0.1,0), "cm"),)

#####################################
load("/Extended/Data/plot_df_sp.RData")
load("/Extended/Data/xi_fr_agg_sp.RData")
load("/Extended/Data/EY_sp.RData")

q <- ggplot(plot_df, aes(x = Time, y = Value, group = Time, color = State)) +
  geom_boxplot(
    size = 0.2, fill = NA, outlier.size = 0.3,
    outlier.shape = 16, outlier.alpha = 0.1
  ) +
  stat_summary(
    fun = mean, geom = "step", aes(group = 1),
    size = 0.1, color = "gray70", linetype = "solid"
  ) +
  geom_point(
    data = data.frame(Time = dates, Value = EY),
    aes(x = Time, y = Value),
    inherit.aes = FALSE, size = 0.2, color = "gray30" #, linetype = "solid"
  ) +
  labs(y = "Expected Strength", title = "Spain") +
  theme_classic() +
  scale_color_gradientn(colors = c("#c20824" , "#ff0080", "#0000ff")) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "italic", size = 24, hjust = 0.5),
    strip.text.x = element_text(size = 16, face = "italic"),
    axis.title = element_text(face = "italic", size = 16),
    axis.title.y = element_text(size = 16, vjust = -0.2),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 16),
    axis.text.x = element_blank(),
    plot.margin = unit(c(0, 0, -2, 0), "cm"),
    axis.line = element_blank(),
    axis.ticks.length.x = unit(0, 'lines'),
    panel.grid = element_blank(),
    legend.title = element_text(face = "italic", size = 16),
    legend.text = element_text(size = 14),
    panel.border = element_rect(colour = "black", fill = NA)
  )
q1<-q

q<-ggplot(xi_fr_agg, aes(x = Time, y = (K+1)-State), col ="red") 
q <- q +  geom_step(size = 0.2, color = "grey70", alpha = 0.8)+ theme_classic()
q <- q +  geom_point(aes(color =  (K+1)-State), size = 0.1, alpha = 0.8)
q<- q+ labs(y = "States") + theme_classic()
q <- q+ scale_color_gradientn(colors=c( "#0000ff", "#ff0080", "#c20824"))
q<- q + theme(legend.position="none", 
              plot.title = element_blank(),
              strip.text.x = element_text(size = 16,face = "italic"),
              axis.title = element_blank(),
              axis.title.y = element_text(size = 16, vjust = -0.2),
              axis.title.x = element_blank(),
              axis.text=element_text(size=5),
              axis.line = element_blank(),
              axis.ticks = element_line(),
              panel.grid = element_blank(),
              #plot.margin = unit(c(-20,0,0,0), "cm"),
              legend.title = element_text(face="italic", size=16),
              legend.text = element_text(size=14),
              panel.border = element_rect(colour = "black", fill = NA))
q2<-q
gv_sp <- q1 / q2 + plot_layout(heights = c(4, 1))
#################

gv_fr<- gv_fr &   theme(axis.title.y = element_text(margin = margin(r = 10)))

gv_it<- gv_it & theme(axis.title.y = element_text(margin = margin(r = 10))) 


hh<- (gv_fr|gv_de)/(gv_it|gv_sp) & theme(legend.position="none", plot.title = element_text(face = "italic",
                                                                                           size = 22),
                                         strip.text.x = element_text(size = 16,face = "italic"),
                                         
                                         axis.title = element_text( face = "italic",size = 18, vjust = -0.8),
                                         axis.text=element_text(size=9, vjust = -0.035),
                                         axis.line = element_blank(),
                                         axis.ticks = element_line(),
                                         panel.grid = element_blank(),
                                         legend.title = element_text(face="italic", size=16),
                                         legend.text = element_text(size=14),
                                         panel.border = element_rect(colour = "black", fill = NA))


hh


ggsave(hh, filename = paste0("Extended/Figures/Figure11.pdf"), units = "cm", width = 16*2, height = 9*2 )
