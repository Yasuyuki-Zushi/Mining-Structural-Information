#####Input your directory#####
#####Input your directory#####
#####Input your directory#####
inputdir <- "C:/.../Code"
codedir <- "C:/.../Code"
#####Input your directory#####
#####Input your directory#####
#####Input your directory#####


#Try Install.packages(), if you don't have them.
library(ggplot2)
library(tsne)

#Input proportion of data for t-sne demo; 0.01 - 1 (all data) 
demo.data.proportion <- 0.01


#color palette 
setwd( inputdir )
color.list <- read.csv("palette.csv")

#loading variants
PubChem_Mw_dataset <- read.csv( "PubChem_Mw.csv"  )
PubChem_Mw_dataset.id <- 1:nrow(PubChem_Mw_dataset)

set.seed(0)
sample.size.tsne <- round( nrow(PubChem_Mw_dataset)*demo.data.proportion )
PubChem_Mw_dataset.id.sample <- sample(PubChem_Mw_dataset.id,size=sample.size.tsne)
PubChem_Mw_dataset.demo <- PubChem_Mw_dataset[PubChem_Mw_dataset.id.sample,]


#t-sne calculation
PubChem_Mw_dataset.demo.tsne <- tsne( PubChem_Mw_dataset.demo[,4:ncol(PubChem_Mw_dataset.demo)],max_iter=1000, k=2, initial_dims=30, perplexity=30 )
PubChem_Mw_dataset.demo.cord <- as.data.frame(PubChem_Mw_dataset.demo.tsne)

#palette setting
exist.superclass.num <- match( unique(PubChem_Mw_dataset.demo$Superclass), color.list$name  )
getPalette <- as.vector( color.list$color[sort(exist.superclass.num)] )

#plot
windows()
ggplot( PubChem_Mw_dataset.demo.cord, mapping = aes(V1,V2,color = PubChem_Mw_dataset.demo$Superclass)) + geom_point() + 
  scale_color_manual(values=getPalette[1:length(unique(PubChem_Mw_dataset.demo$Superclass))],name="Superclass" )

