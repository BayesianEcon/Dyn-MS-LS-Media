

rm(list = ls())

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

#install all the packages necessary to run the scripts
# Compatibility tested with R 2020-10-10, "Bunny-Wunnies Freak Out"

#For Windows users: if source compilation fails try install the package individually with:
#install.packages("packagename", type = "binary")
#If you run rcpp scripts from R makes sure that Rtools is correctly installed
#https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html

#For Mac users: if you experience troubles in the installation of Rcpp with an M1 processor 
#(and later versions), please refer to the following link for further instructions
#https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/

#packages for running the models

list.of.packages <- c("mvtnorm", "Rcpp", "RcppDist", "RcppParallel", "RcppArmadillo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load libraries to check that everything is running smoothly

library(mvtnorm)
library(Rcpp)
library(RcppDist)
library(RcppParallel)
library(RcppArmadillo)

#packages for data manipulation and plotting part 1

list.of.packages <- c( "ggplot2", "ggrepel", "tidyr", "scales", "dplyr", "data.table", "patchwork", "ggpubr", "metR", "igraph")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")

#load libraries to check that everything is running smoothly

library(ggplot2)
library(ggrepel)
library(tidyr)
library(scales)
library(dplyr)
library(data.table)
library(patchwork)
library(ggpubr)
library(metR)
library(igraph)

#packages for data manipulation and plotting part 2

list.of.packages <- c("lubridate","grid", "ggraph","ggiraphExtra","ggtext","ggnewscale","tidygraph", "ggpmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, type = "binary")

#load libraries to check that everything is running smoothly

library(lubridate)
library(grid)
library(ggraph)
library(ggiraphExtra)
library(ggtext)
library(ggnewscale)
library(tidygraph)
library(ggpmisc)

#packages for data manipulation and plotting part 2

list.of.packages <- c("udpipe", "quanteda", "quanteda.textstats", "plm", "readtext","readr", "reshape2","stringi","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load libraries to check that everything is running smoothly

library(udpipe)
library(quanteda)
library(quanteda.textstats)
library(plm)
library(readtext)
library(readr)
library(plyr)
library(reshape2)
library(stringi)
library(stringr)


rm(list = ls())

#packages for data manipulation and plotting part 3

list.of.packages <- c("MCMCpack", "tidyverse", "gtable")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#load libraries to check that everything is running smoothly

library(MCMCpack)
library(tidyverse)
library(gtable)
