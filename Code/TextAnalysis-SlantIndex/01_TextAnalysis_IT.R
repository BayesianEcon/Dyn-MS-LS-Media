rm(list = ls())

#load the libraries

list.of.packages <- c("udpipe", "quanteda", "quanteda.textstats", "plm", "readtext","readr",
                      "plyr","reshape2", "lubridate","ggplot2","ggrepel","stringi","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

########CHANGE YOUR PATH ###########
setwd("~/Desktop/Repository/")
#####################################

library(udpipe)
library(quanteda)
library(quanteda.textstats)
library(plm)
library(readtext)
library(readr)
library(plyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(stringi)
library(stringr)



quanteda_options(verbose = TRUE)
quanteda_options(threads = 4)

load("Data/TextAnalysis-SlantIndex/Italy/outlets1.RData")
load("Data/TextAnalysis-SlantIndex/Italy/outlets2.RData")
load("Data/TextAnalysis-SlantIndex/Italy/outlets3.RData")
load("Data/TextAnalysis-SlantIndex/Italy/outlets4.RData")

outlets<-cbind(outlets1,outlets2,outlets3,outlets4)
rm(outlets1,outlets2,outlets3,outlets4)
gc()

source("Code/TextAnalysis-SlantIndex/stopmore.R")

ng<- 1:2

party_names<- c("Centro Democratico.txt", "Forza Italia.txt", "Fratelli d'Italia.txt", "Lega - Salvini Premier.txt", "MoVimento 5 Stelle.txt", 
                "Partito Democratico.txt", "Scelta Civica.txt", "Sinistra Ecologia Libertà - Milano.txt" , "Udc Italia.txt" )

orientation<-data.frame(Var2 = party_names, score = c( 5 , 6.71 , 7.86 ,8.86 , 4.67 , 3.57, 5.42, 1.28, 5.28 )) 		


#create pages scripts

dom <- function(x) strsplit(gsub("http://|https://|www\\.", "", x), "/")[[c(1, 1)]]

party_posts <- read_csv("Data/TextAnalysis-SlantIndex/Italy/party_posts.csv")

party_posts<-party_posts[party_posts$`Page Name` != "SVP - Südtiroler Volkspartei",]
party_posts$domain <- sapply(party_posts$Link, FUN= dom)
party_posts$domain <- gsub(".com|.it|.org", "",party_posts$domain )
party_posts$text<- paste(party_posts$Message,party_posts$`Link Text`, party_posts$Description)
party_posts$date <-   as_date(party_posts$Created)


outlets$date<-as.numeric(as_date(outlets$date))
outlets<-outlets[outlets$date >= 16436 &  outlets$date <= 17166  ,]

outlets$text<- paste(outlets$Message, outlets$`Link Text`, outlets$Description)
outlets$text <- gsub("@", "", outlets$text)
outlets$text <- gsub("https?://.+", "", outlets$text)
outlets$text <- gsub("\\d+\\w*\\d*", "", outlets$text)

outlets$date <-   as_date(outlets$Created)

# CorpusOu <-corpus(outlets$text)
# 
# mytokenso <- tokens(CorpusOu, remove_url = TRUE, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_separators = TRUE,  what = "word", verbose = T )
# mytokenso <- tokens_tolower(mytokenso)
# mytokenso <- tokens_select(mytokenso, pattern = c(stopmore, stopwords(language = "it"), stopwords("en") ), selection = "remove")
# 
# # #toks_ngramo <- tokens_ngrams(mytokenso, n = ng)
# # 
# # rm(mytokenso)
# # rm(toks_ngramo)
# dfm_o <- dfm(mytokenso)
# 
#list <- featnames(dfm_o)
# 
#save(list, file = "dictionary_outlets.RData")
# rm(dfm_o)

load("Data/TextAnalysis-SlantIndex/Italy/dictionary_outlets.RData")


party_posts$text<- gsub("Passa alla slide successiva", "",party_posts$text)


party_posts$text <- gsub("@", "", party_posts$text)
party_posts$text <- gsub("https?://.+", "", party_posts$text)
party_posts$text  <- gsub("\\d+\\w*\\d*", "", party_posts$text)
party_posts$text <- gsub("[[:punct:]]", " ", party_posts$text )
party_posts$text <- gsub("Timeline Photos NA", "", party_posts$text )


# Remove spaces and newlines
party_posts$text  <- gsub("\n", " ", party_posts$text )
party_posts$text <- gsub("^\\s+", "", party_posts$text )
party_posts$text <- gsub("\\s+$", "", party_posts$text )
party_posts$text  <- gsub("[ |\t]+", " ",party_posts$text )
party_posts$text  <- gsub("\\.", " ",party_posts$text )


dates<-unique(outlets$date)
dates<- dates[order(dates)]

V<-data.frame(mean= 0, sd = 0, nseed = 0, date = 0)

path1<- "Data/TextAnalysis-SlantIndex/Italy/party/"
path2<- "Data/TextAnalysis-SlantIndex/Italy/outlet/"

party_posts_sub<-party_posts

AGG<-aggregate(party_posts_sub$text, by = list(party_posts_sub$`Page Name`), FUN = paste) 

for(i in 1:nrow(AGG)){
  write.table(AGG[i,2], file = paste0("Data/TextAnalysis-SlantIndex/Italy/party/", AGG[i,1], ".txt"), row.names = F)
  #print(i)
}

Pol_txt<- readtext("Data/TextAnalysis-SlantIndex/Italy/party/*.txt",  cache = FALSE)

Corpus <-corpus(Pol_txt)

mytokens <- tokens(Corpus, remove_url = TRUE, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_separators = TRUE,  what = "word", verbose = T )
mytokens <- tokens_tolower(mytokens)
mytokens <- tokens_select(mytokens, pattern = list, selection = "keep")
mytokens <- tokens_select(mytokens, pattern = c(stopmore, stopwords(language = "it"), stopwords("en") ), selection = "remove")

toks_ngram <- tokens_ngrams(mytokens, n = ng)

mydfm <- dfm(toks_ngram, keep = list, remove = c( stopwords("en"), stopwords(language = "it"), stopmore))
mydfm <- dfm_select(mydfm, pattern = list)

mydfm<- dfm_keep(mydfm, min_nchar = 3)


tfidf<-dfm_tfidf(mydfm)

keyw<- 100 #ottimale

if( "Centro Democratico.txt" %in% rownames(tfidf)){
  CD<-convert(t(tfidf[rownames(tfidf) == "Centro Democratico.txt"  ,]),to = "data.frame")
  CD<-CD[order(CD$`Centro Democratico.txt`, decreasing=T),]
  CD<-CD[CD$`Centro Democratico.txt`> 0, ]
  #CD<-CD[1:(round(keyw*nrow(CD))),]
  CD<-CD[1:keyw,]
  
  
  }else{CD <- data.frame() ; CD$doc_id <-character() }

if( "Forza Italia.txt" %in% rownames(tfidf)){
  PDL<-convert(t(tfidf[rownames(tfidf) == "Forza Italia.txt",]),to = "data.frame")
  PDL<-PDL[order(PDL$`Forza Italia.txt`, decreasing=T),]
  PDL<-PDL[PDL$`Forza Italia.txt`> 0, ]
  #PDL<-PDL[1:(round(keyw*nrow(PDL))),]
  PDL<-PDL[1:keyw,]
  
  }else{PDL <- data.frame() ; PDL$doc_id <-character() }

if( "Fratelli d'Italia.txt" %in% rownames(tfidf)){
  FDI<-convert(t(tfidf[rownames(tfidf) == "Fratelli d'Italia.txt",]),to = "data.frame")
  FDI<-FDI[order(FDI$`Fratelli d'Italia.txt`, decreasing=T),]
  FDI<-FDI[FDI$`Fratelli d'Italia.txt`> 0, ]
  #FDI<-FDI[1:(round(keyw*nrow(FDI))),]
  FDI<-FDI[1:keyw,]
  
  }else{FDI <- data.frame() ; FDI$doc_id <-character() }

if( "Lega - Salvini Premier.txt" %in% rownames(tfidf)){
  LEGA<-convert(t(tfidf[rownames(tfidf) == "Lega - Salvini Premier.txt" ,]),to = "data.frame")
  LEGA<-LEGA[order(LEGA$`Lega - Salvini Premier.txt`, decreasing=T),]
  LEGA<-LEGA[LEGA$`Lega - Salvini Premier.txt`> 0, ]
  #LEGA<-LEGA[1:(round(keyw*nrow(LEGA))),]
  LEGA<-LEGA[1:keyw,]
  
  }else{LEGA <- data.frame() ; LEGA$doc_id <-character() }

if( "MoVimento 5 Stelle.txt" %in% rownames(tfidf)){
  M5S<-convert(t(tfidf[rownames(tfidf) == "MoVimento 5 Stelle.txt",]),to = "data.frame")
  M5S<-M5S[order(M5S$`MoVimento 5 Stelle.txt`, decreasing=T),]
  M5S<-M5S[M5S$`MoVimento 5 Stelle.txt`> 0, ]
  #M5S<-M5S[1:(round(keyw*nrow(M5S))),]
  M5S<-M5S[1:keyw,]
  
  
  }else{M5S <- data.frame() ; M5S$doc_id <-character() }

if( "Partito Democratico.txt" %in% rownames(tfidf)){
  PD<-convert(t(tfidf[rownames(tfidf) == "Partito Democratico.txt" ,]),to = "data.frame")
  PD<-PD[order(PD$`Partito Democratico.txt`, decreasing=T),]
  PD<-PD[PD$`Partito Democratico.txt`> 0, ]
  #PD<-PD[1:(round(keyw*nrow(PD))),]
  PD<-PD[1:keyw,]
  
  }else{PD <- data.frame() ; PD$doc_id <-character() }

if( "Scelta Civica.txt" %in% rownames(tfidf)){
  SC<-convert(t(tfidf[rownames(tfidf) == "Scelta Civica.txt",]),to = "data.frame")
  SC<-SC[order(SC$`Scelta Civica.txt`, decreasing=T),]
  SC<-SC[SC$`Scelta Civica.txt`> 0, ]
  #SC<-SC[1:(round(keyw*nrow(SC))),]
  SC<-SC[1:keyw,]
  
  }else{SC <- data.frame() ; SC$doc_id <-character() }

# SVP<-convert(t(tfidf[60,]),to = "data.frame")
# SVP<-SVP[order(SVP$`SVP - Südtiroler Volkspartei .txt`, decreasing=T),]
# SVP<-SVP[1:82,]

if( "Udc Italia.txt" %in% rownames(tfidf)){
  UDC<-convert(t(tfidf[rownames(tfidf) == "Udc Italia.txt" ,]),to = "data.frame")
  UDC<-UDC[order(UDC$`Udc Italia.txt`, decreasing=T),]
  UDC<-UDC[UDC$`Udc Italia.txt`> 0, ]
  #UDC<-UDC[1:(round(keyw*nrow(UDC))),]
  UDC<-UDC[1:keyw,]
  
  }else{UDC <- data.frame() ; UDC$doc_id <-character() }

if( "Sinistra Ecologia Libertà - Milano.txt" %in% rownames(tfidf)){
  VENDOLA<-convert(t(tfidf[rownames(tfidf) == "Sinistra Ecologia Libertà - Milano.txt" ,]),to = "data.frame")
  VENDOLA<-VENDOLA[order(VENDOLA$`Sinistra Ecologia Libertà - Milano.txt`, decreasing=T),]
  VENDOLA<-VENDOLA[VENDOLA$`Sinistra Ecologia Libertà - Milano.txt`> 0, ]
  #VENDOLA<-VENDOLA[1:(round(keyw*nrow(VENDOLA))),]
  VENDOLA<-VENDOLA[1:keyw,]
  
  }else{VENDOLA <- data.frame() ; VENDOLA$doc_id <-character() }


key_features<-c(CD$doc_id,
                PDL$doc_id,
                FDI$doc_id,
                LEGA$doc_id,
                M5S$doc_id,
                PD$doc_id,
                SC$doc_id,
                UDC$doc_id,
                VENDOLA$doc_id
)


key_features<-key_features[!is.na(key_features)]

dict<-dictionary(list(key= key_features))

mydfm<-dfm_select(mydfm, pattern = dict, case_insensitive = FALSE, selection = "keep")

DBplane<-data.frame( fb_name=0,ncomments= 0, leaning = 0, date =0 )

colnames(outlets)[7]<- "comments_count"
colnames(outlets)[1]<- "from_name"

outlet_list<- unique(outlets$from_name)

dates<- unique(outlets$date)

for(k in 1:length(dates)){
  
  # 
  junk2 <- dir(path= path2,  pattern="*.txt")
  # 
  if(length(junk2)!= 0) { file.remove(paste0(path2, junk2))}
  
  
  outlets_sub<-outlets[outlets$date == dates[k],]
  
  outlets_comments_db<-aggregate(outlets_sub$comments_count, by = list(outlets_sub$from_name), FUN = sum)
  colnames(outlets_comments_db)<-c("fb_name", "ncomments")
  
  missing<-outlet_list[!(outlet_list %in%   outlets_comments_db$fb_name)]
  if(length(missing)!= 0){
    block_missing<- data.frame(fb_name =missing, ncomments = 0 )
    outlets_comments_db<-rbind(outlets_comments_db, block_missing )
    rm(missing, block_missing)
    }
  
  AGG_out<-aggregate(outlets_sub$text, by = list(outlets_sub$from_name), FUN = paste)
  

  for(i in 1:nrow(AGG_out)){
    write.table(AGG_out[i,2], file = paste0("Data/TextAnalysis-SlantIndex/Italy/outlet/",AGG_out[i,1], ".txt"), row.names = F)
    #print(i)
  }
  
  #create corpus
  
  TXT<- readtext("Data/TextAnalysis-SlantIndex/Italy/outlet/*.txt",  cache = FALSE)
  
  Corpus_outlet <-corpus(TXT)
  
  tokens_outlet <- tokens(Corpus_outlet, remove_url = TRUE, remove_punct = TRUE, remove_numbers = TRUE, remove_symbols = TRUE, remove_separators = TRUE,  what = "word", verbose = T )
  tokens_outlet <- tokens_tolower(tokens_outlet)
  tokens_outlet <- tokens_select(tokens_outlet, pattern = c(stopmore, stopwords(language = "it"), stopwords("en") ), selection = "remove")
  toks_ngram_outlet <- tokens_ngrams(tokens_outlet, n = ng)
  
  dfm <- dfm(toks_ngram_outlet,tolower = FALSE)
  dfm <- dfm_select( dfm, pattern = dict,  selection = "keep",  case_insensitive = TRUE  )
  # dfm<-dfm_match(dfm , features = featnames(mydfm) )
  mydfm_s<- dfm_match(mydfm , features = featnames(dfm) )
  
  tstat <- textstat_simil(dfm, mydfm_s , method = "cosine", margin = "documents")

  M<-as.matrix(tstat)
  M[is.na(M)]<-0
  
  DB<-reshape2::melt(M)
  DB<-na.omit(DB)
  DB$Var2<-as.character(DB$Var2)
  DB<-join(DB,orientation, by = "Var2" )

  colnames(DB) <- c("page", "party", "value", "score")
  
  DB$ones<- rnorm(nrow(DB))
 
  mod1<- lm(formula = value ~ page + party, data = DB )
  summary(mod1)
  
  DB$residuals<- mod1$residuals
  DB$score_res<-DB$score*DB$residuals
  DB<-na.omit(DB)
  
  spectrum<- aggregate(DB$score_res, by = list(DB$page), FUN = sum)
  
  colnames(spectrum)<-c("fb_name", "leaning")

  
  spectrum$fb_name<-gsub(".txt","",spectrum$fb_name)
  
  ###### add missing
  missing<-outlet_list[!(outlet_list %in%   spectrum$fb_name)]
  if(length(missing)!= 0){
  block_missing<- data.frame(fb_name =missing, leaning = 0 )
  spectrum<-rbind(spectrum, block_missing )}
  #####
  
  outlets_comments_db<-join(outlets_comments_db, spectrum , by = "fb_name")
  outlets_comments_db$date<-dates[k]
  
  DBplane<-rbind(DBplane, outlets_comments_db)
  
  print(k)
}


DBplane<-DBplane[-1,]

colnames(DBplane)[1]<-"i"
colnames(DBplane)[4] <- "t"

DBplane$i<-gsub(" ", ".", DBplane$i)
DBplane$i<-gsub("'", ".", DBplane$i)
DBplane$i<-gsub("\\.-\\.", ".", DBplane$i)


save(DBplane, file = "Data/TextAnalysis-SlantIndex/Italy/DBplane_check_it.RData")

load("Data/TextAnalysis-SlantIndex/Italy/DBplane_check_it.RData")


DBplane_new<- DBplane

colnames(DBplane_new)[1]<-"fb_name"
colnames(DBplane_new)[2] <- "ncomments"
colnames(DBplane_new)[3] <- "leaning"
colnames(DBplane_new)[4] <- "t"


agg<- aggregate(DBplane_new$leaning, by = list(DBplane_new$fb_name), FUN = mean)
agg<-agg[order(agg$x),]
agg$Group.1

load("Data/Dynamic/DataEnv_IT_all.RData")

page_list<- unique(DBplane$fb_name)
page_list[13]<- "Il.Messaggero"

DBplane_new$i<- as.numeric(as.factor(DBplane_new$fb_name))
table(DBplane_new$fb_name)
table(DBplane_new$t)


rownames(DBplane_new)<-1:nrow(DBplane_new)
rownames(DBplane)<-1:nrow(DBplane)



DBplane_new<- DBplane_new[DBplane_new$fb_name %in% page_list, ]
DBplane_new<-DBplane_new[DBplane_new$t %in% unique(DBplane$t), ]

agg<- aggregate(DBplane$leaning, by = list(DBplane$fb_name), FUN = mean)
agg<-agg[order(agg$x),]
agg$Group.1

DBplane$fb_name<-factor(DBplane$fb_name, levels = agg$Group.1)
DBplane$leaning<- DBplane$leaning - mean(DBplane$leaning )


p<- ggplot(DBplane, aes(x=leaning, y= fb_name, group = fb_name )) 
p<- p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),   geom="errorbar", color="red", width=0.2) 
p<- p + stat_summary(fun=mean, geom="point", color="red", size = 3) 
p<- p + labs(x = "Text-Analysis Leaning", y = "" )
p<- p +  theme_minimal()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                    size = 36, hjust = 0.5),
                                strip.text.x = element_text(size = 24,face = "bold"),
                                
                                axis.title = element_text( face = "bold",size = rel(1)),
                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
                                axis.title.x = element_text(size = 22, vjust = -0.2),
                                axis.text=element_text(size=16),
                                axis.line = element_blank(),
                                axis.ticks = element_line(),
                                legend.title = element_text(face="italic", size=16),
                                legend.text = element_text(size=14),
                                panel.border = element_rect(colour = "black", fill = NA))


p_it<-p

ggsave(p, file = "Figures/TextAnalysis-SlantIndex/OutletsLeaning_IT.pdf", width = 16, height = 9) 


party_names<- c("Centro Democratico.txt", "Forza Italia.txt", "Fratelli d'Italia.txt", "Lega - Salvini Premier.txt", "MoVimento 5 Stelle.txt", 
                "Partito Democratico.txt", "Scelta Civica.txt", "Sinistra Ecologia Libertà - Milano.txt" , "Udc Italia.txt" )



pn <- c(
  "CD",
  "PDL",
  "FDI",
  "LEGA",
  "M5S",
  "PD",
  "SC",
  "SEL",
  "UDC"
)

pl<- orientation$score
DB_party<- data.frame(party_name = pn, party_leaning = pl)

q<-ggplot(DB_party, aes(x = party_leaning , y = 0))+ geom_hline(yintercept = 0)+geom_point(color = "red", size = 7) +geom_text_repel(aes(label = party_name), direction = "x",nudge_y = 0.09,  size = 7, min.segment.length = Inf) + ylim(-0.5,0.5)
q<- q + labs(x = "Party Leaning", y = "" ) + xlim(1,11)
q<- q +  theme_minimal()+ theme(legend.position="bottom", plot.title = element_text(face = "bold",
                                                                                    size = 36, hjust = 0.5),
                                strip.text.x = element_text(size = 24,face = "bold"),
                                
                                axis.title = element_text( face = "bold",size = rel(1)),
                                axis.title.y = element_text(size = 22, angle=90,vjust =2),
                                axis.title.x = element_text(size = 22, vjust = -0.2),
                                axis.text=element_text(size=16),
                                axis.text.y=element_blank(),
                                axis.line = element_blank(),
                                axis.ticks = element_line(),
                                axis.ticks.y = element_blank(),
                                panel.grid = element_blank(),
                                legend.title = element_text(face="italic", size=16),
                                legend.text = element_text(size=14),
                                panel.border = element_rect(colour = "black", fill = NA))


q_it<-q




ggsave(q, file = "Figures/TextAnalysis-SlantIndex/PartiesLeaning_IT.pdf", width = 16, height = 4.5)
