#####Input your directory#####
#####Input your directory#####
#####Input your directory#####
inputdir <- "C:/.../Input"
codedir <- "C:/.../Code"
#####Input your directory#####
#####Input your directory#####
#####Input your directory#####


#Try Install.packages(), if you don't have them.
library(caret)
library(xgboost)
library(Matrix)



#loading variants
setwd(inputdir)
PubChem_Mw_dataset <- read.csv( "PubChem_Mw.csv"  )

na.col <- which( is.na(colSums(PubChem_Mw_dataset[,4:ncol(PubChem_Mw_dataset)])) )
if( sum(na.col)==0 ){
  PubChem_Mw_dataset_prc <- PubChem_Mw_dataset
} else {
  rem.col <- match( names(na.col), colnames(PubChem_Mw_dataset) )
  PubChem_Mw_dataset_prc <- PubChem_Mw_dataset[,-rem.col]
}
all.data.final <- PubChem_Mw_dataset_prc[,4:ncol(PubChem_Mw_dataset_prc)]

#label data
Mwdata <-  log10( PubChem_Mw_dataset_prc[,3] )
Mwdata.id <- 1:length(Mwdata)


#data split for test data
set.seed(0)
sample.size.test <- round( length(Mwdata)*0.1 )
index.test <- Mwdata.test.id <- sample(Mwdata.id,size=sample.size.test)

no.data.test <- which( index.test > length(Mwdata) )
if( sum(no.data.test)==0){
  index.test <- index.test
} else {
  index.test <- index.test[-no.data.test]
}

all.data.final.test <- as.matrix( all.data.final[index.test,] )
labels.test <- as.numeric(Mwdata[index.test])


#data split for validation data
sample.size.validation <- round( length(Mwdata)*0.1 )
index.validation <- Mwdata.validation.id <- sample(Mwdata.id[-index.test],size=sample.size.validation)

no.data.validation <- which( index.validation > length(Mwdata) )
if( sum(no.data.validation)==0 ){
  index.validation <- index.validation
} else {
  index.validation <- index.validation[-no.data.validation]
}

all.data.final.validation <- as.matrix( all.data.final[index.validation,] )
labels.validation <- as.numeric(Mwdata[index.validation])



#data split for train data
index.train <- Mwdata.id[-c(index.test,index.validation)]
all.data.final.train <- as.matrix( all.data.final[index.train,] )
labels.train <- as.numeric( Mwdata[index.train] )


###data convert for xgb
data.test <- all.data.final.test
data.test[which(is.na(data.test))] <- 0

data.validation <- all.data.final.validation
data.validation[which(is.na(data.validation))] <- 0

data.train <- all.data.final.train
data.train[which(is.na(data.train))] <- 0

bind.data.train <- as.data.frame( cbind(labels.train,data.train) )
bind.data.test <- as.data.frame( cbind(labels.test,data.test) )
bind.data.validation <- as.data.frame( cbind(labels.validation,data.validation) )

colnames(bind.data.train) <- c( "labels.train",paste("fp",1:(ncol(data.train)),sep="") )
colnames(bind.data.test) <- c( "labels.test",paste("fp",1:(ncol(data.test)),sep="") )
colnames(bind.data.validation) <- c( "labels.validation",paste("fp",1:(ncol(data.validation)),sep="") )

options(na.action='na.pass')
train_dmat <- xgb.DMatrix(
  sparse.model.matrix(labels.train~.-1, data = bind.data.train),
  label = bind.data.train$labels.train
)
options(na.action='na.omit')

options(na.action='na.pass')
test_dmat <- xgb.DMatrix(
  sparse.model.matrix(labels.test~.-1, data = bind.data.test),
  label = bind.data.test$labels.test
)
options(na.action='na.omit')

options(na.action='na.pass')
validation_dmat <- xgb.DMatrix(
  sparse.model.matrix(labels.validation~.-1, data = bind.data.validation),
  label = bind.data.validation$labels.validation
)
options(na.action='na.omit')



#xgb paramter
l_params = list(booster = 'gbtree', 
                objective = 'reg:linear',
                eval_metric = 'rmse',
                eta = 0.02,
                max_depth = 7,
                min_child_weight = 1,
                alpha = 1,
                lambda = 1
) 

#xgb Dmatrix
xgb_model <- xgb.train(
  data = train_dmat,
  nrounds = 10000,
  params = l_params,
  watchlist = list(train = train_dmat),
  early_stopping_rounds = 100
)


#xgb calculation
xgb_pred <- predict(xgb_model,validation_dmat)


#RMSE
rmse.calc <- round( sqrt( sum((bind.data.validation$labels.validation - xgb_pred)^2) / length(xgb_pred) ),digits=3)


#plot range
wide.fa <- 1.5
tiny.v <- .Machine$double.eps
if( min(bind.data.validation$labels.validation)<0 ){ xlim.min <- min(bind.data.validation$labels.validation) * wide.fa } else 
{ xlim.min <- min(bind.data.validation$labels.validation) - wide.fa*(min(bind.data.validation$labels.validation) + tiny.v) }
if( max(bind.data.validation$labels.validation)<0 ){ xlim.max <- max(bind.data.validation$labels.validation) - wide.fa*max(bind.data.validation$labels.validation) } else 
{ xlim.max <- (max(bind.data.validation$labels.validation) + tiny.v) * wide.fa }

if ( min(xgb_pred)<0 ){ ylim.min <- min(xgb_pred) * wide.fa } else
{ ylim.min <- min(xgb_pred) - wide.fa * (min(xgb_pred) + tiny.v) }
if ( max(xgb_pred)<0 ){ ylim.max <- max(xgb_pred) - wide.fa * max(xgb_pred) } else
{ ylim.max <- (max(xgb_pred) + tiny.v) * wide.fa }


#plot
windows(width = 100, height = 100)
plot(bind.data.validation$labels.validation, xgb_pred,xlab="Data for validation",ylab="Prediction",xlim=c(xlim.min, xlim.max),ylim=c(ylim.min,ylim.max) )
lm.res <- lm(xgb_pred ~ bind.data.validation$labels.validation)
intercep.v <- round(summary(lm.res)[4]$coefficients[1],digits=3)
slop.v <- round(summary(lm.res)[4]$coefficients[2],digits=3)
r2.v <- round(as.numeric(summary(lm.res)[8]),digits=3)
abline(a=intercep.v,b=slop.v,col="red")
abline(a=0, b=1, col="blue")
abline(a=1, b=1, col="blue",lty=2)
abline(a=-1, b=1, col="blue",lty=2)
legend("topleft",c( paste("Slope=",slop.v,", Intercept=",intercep.v,", R-sqr=",r2.v,sep=""),
                    paste("RMSE=",round(rmse.calc,digits=2),", n=",length(xgb_pred),sep="") ) )

