library(matrixStats)
library(dplyr)
library(timeSeries)
library(StatMatch)
library(FNN)
library(diffusionMap)
library(rARPACK)
library(leaps)
library(energy)
library(pracma)
library(rowr)


currencies <- read.csv(file = "currencies.csv", header = TRUE)
commodities <- read.csv(file = "commodities.csv", header = TRUE)
complete <- merge(currencies, commodities, by='Date', sort = FALSE)
complete <- complete[complete.cases(complete),]
currencies <- complete[,2:ncol(currencies)]
commodities <- complete[,(ncol(currencies)+3):ncol(complete)]


##################################

## PCA

##################################

run_PCA = function(X_train, nrFactors){
  
  myPCA = prcomp(X_train, center = TRUE, scale. = TRUE, retx = TRUE);
  pcaRot = myPCA$rotation
  EIGS_PCA = myPCA$x
  
  q = myPCA$sdev
  
  FACTORS = EIGS_PCA[, 1:nrFactors];
  
  
  return(FACTORS);
}

#########################

## PLOT FUNCTION

#########################

plot_cumsum = function(CUM_SUM_MTX,MAIN,low,high,location){
  # INPUT: CUM_SUM_MTX, a matrix of size nrDays  x nrTickers 
  #   which plots the time series (one from each column) on the same plot.
  
  nrTickers = dim(CUM_SUM_MTX)[2];
  nrDays = dim(CUM_SUM_MTX)[1];
  myColors = c('black','red','blue', 'green','yellow','magenta')
  
  # pdf(file = 'cumulative_pnl_across_time.pdf', height = 12,  width = 8);
  # plot( CUM_SUM_MTX[ ,1], col = myColors[1], type ='l', pch=3, ylim = c(-2,6), lwd=2);
  plot (c(0, nrDays),c(low,high),type='n', xlab ='days', ylab = 'Prices',main =  MAIN) # sets the x and y axes scales
  for ( i in 1 : nrTickers){
    # plots the time series on the same graph
    lines( 1: nrDays, CUM_SUM_MTX[ ,i], col = myColors[i], lwd=2);  # , ylim = c(-1,1)
  }
  
  legend( location, legend = colnames(CUM_SUM_MTX), lty = 1, lwd = 2, col = myColors)
}




nrDays<-nrow(complete)
table<-NULL
windowdata<- complete[,-1] #getting rid of date column 
complete_dropdates<-windowdata
estimated_values <- data.frame(matrix(nrow = (nrDays-101), ncol = 0))

###########

## LINEAR REGRESSION BY CORRELATIONS

###########

for(j in 1:38){
  
  estimate <- data.frame()

  for(d in 101:(nrDays-1)){
    
    response <- complete_dropdates[,j]
    predict <- complete_dropdates[,]
    response <- response[(d-100):(d-1)]
    predict <- predict[(d-99):d,]
    
    correlation <- cor(as.matrix(complete_dropdates), method = "pearson")
    diag(correlation) <- 0
    abs_correlation <- abs(correlation)
    
    y <- sort(abs_correlation[,j], decreasing = TRUE)[1:10]
    pos <- which(abs_correlation[,j] %in% y)
    model.elements <- predict[,pos]
    
    model.data <- data.frame(response,predict)
    model <- lm(response~as.matrix(model.elements))
    v <- model$coefficients
    est <- as.matrix(complete_dropdates[d,pos]) %*% v[2:11] + v[1]
    estimate <- rbind(estimate,est)
    
  }
  
estimated_values<-data.frame(estimated_values,estimate)
  
}
estimated_values <- as.matrix(estimated_values)

#####

## LINEAR REGRESSION BY STEPWISE METHOD

## NOTE: NOT CHANGED FROM 12/12/2015

####

estimated_values_stepwise <- data.frame(matrix(nrow = (nrDays-101), ncol = 0))

for(j in 1:38){
  
  estimate <- data.frame()
  
  for(d in 101:(nrDays-1)){
    
    response <- complete_dropdates[,j]
    predict <- complete_dropdates[,]
    response <- response[(d-100):(d-1)]
    predict <- predict[(d-99):(d),]
    
    #compute correlations for the relevant time period 
    model.data <- data.frame(response, predict)
    backwards <- regsubsets(response~.,data=model.data,method= "forward", nvmax = 40)
    topten <- summary(backwards)$which[10,]
    model.elements <- predict[,topten[2:39]]
    
    model<-lm(response~as.matrix(model.elements))
    v <- model$coefficients
    est <- as.matrix(complete_dropdates[d,topten[2:39]]) %*% v[2:11] + v[1]
    estimate <- rbind(estimate,est)
  }
  
  estimated_values_stepwise <- data.frame(estimated_values_stepwise,estimate)
  
}



###################################

######## PCA #########

###################################

nrLags = 1
lookback = 100; # sliding window size in days
nrTickers <- ncol(windowdata)
knn_Pred<- NULL; lm_Pred <- NULL;
knn_Pred = array(0, dim = c(nrDays,nrTickers)); # Creating array of predicted values using K nearest neighbor
lm_Pred = array(0, dim = c(nrDays,nrTickers)); # Creating array of predicted values using linear regression

for(i in (lookback + nrLags) : (nrDays - nrLags)) {
  xtrain <- windowdata[ (i-lookback) : (i-1), ];
  xtest <- windowdata[ i, , drop=FALSE]
  
  both <- rbind(xtrain , xtest)
  
  ans <-run_PCA(both, 10) # perform PCA on the 101xnrFeatures matrix

  xtrain<-ans[1:lookback,]
  xtest<-ans[ (lookback + 1), ]
  
  for(j in 1:38){
    y_train = windowdata [ (i-lookback+1) : i,j]
    y_test = windowdata [i+1,j ]
    
    # k Nearest Neighbor section and entry into matrix
    
    y_hat = knn.reg(train = xtrain, test = xtest, y = y_train, k = 10);
    y_hat = y_hat$pred;
    
    knn_Pred[i,j] = y_hat
    
    # Linear Regression section and entry into matrix
    
    model<-lm(y_train~xtrain)
    v <- model$coefficients
    est <- xtest %*% v[2:11] + v[1]
    
    lm_Pred[i,j] = est
    
  }
}


#Get rid of zeroes at beginning of matrix
knn_Pred = knn_Pred[(lookback+1):(nrDays-1),]
lm_Pred = lm_Pred[(lookback+1):(nrDays-1),]

# Use the "which" function to find specific indices of columns with currency names

USD_pos <- which (colnames(windowdata)=="USD")
BRL_pos <- which (colnames(windowdata)=="BRL")
AUD_pos <- which (colnames(windowdata)=="AUD")
RUB_pos <- which (colnames(windowdata)=="RUB")
ARS_pos <- which (colnames(windowdata)=="ARS")
Oil_pos <- which (colnames(windowdata)=="Oil")


# Vector that has 6 numbers that, in the loop, we use index of
posvector <-c(USD_pos,BRL_pos,AUD_pos,RUB_pos,ARS_pos,Oil_pos)
errormatrix <- matrix(nrow = 4, ncol = 6) 
rownames(errormatrix) <- c("knnPCA","lmPCA","corLM","stepLM")
colnames(errormatrix) <- c("USD","BRL","AUD","RUB","ARS","Oil")

sumerrormatrix <- matrix(nrow= 4, ncol = 6)
rownames(sumerrormatrix) <- c("knnPCA","lmPCA","corLM","stepLM")
colnames(sumerrormatrix) <- c("USD","BRL","AUD","RUB","ARS","Oil")

predictionmatrix <- matrix(nrow = 1084,ncol = 30)

par(mfrow=c(3,2))

j <- 1
for(i in 1:length(posvector)){
  currentPos<-posvector[i] # get index from vector of positions we want
  knnPCA <- knn_Pred[,currentPos]
  lmPCA <- lm_Pred[,currentPos]
  corLM <- estimated_values[,currentPos]
  stepLM <- estimated_values_stepwise[,currentPos]
  realCurrency <- windowdata[(102:1185),currentPos]
  
  
  predictionmatrix[,j] <- knnPCA; 
  predictionmatrix[,j+1] <- lmPCA; 
  predictionmatrix[,j+2] <- corLM; 
  predictionmatrix[,j+3] <- stepLM; 
  predictionmatrix[,j+4] <- realCurrency; 
  


  # Compute Errors and put into matrix
  
  Error_knnPCA <- (realCurrency - knnPCA)^2
  Error_knnPCA <- sqrt(Error_knnPCA)
  Error_lmPCA <- (realCurrency - lmPCA)^2
  Error_lmPCA <- sqrt(Error_lmPCA)
  Error_corLM <- (realCurrency - corLM)^2
  Error_corLM <- sqrt(Error_corLM)
  Error_stepLM <- (realCurrency - stepLM)^2
  Error_stepLM <- sqrt(Error_stepLM)
  
  knnPCA_mean <- mean(Error_knnPCA); errormatrix[1,i] <- knnPCA_mean; sumerrormatrix[1,i] <- sum(Error_knnPCA)
  lmPCA_mean <- mean(Error_lmPCA); errormatrix[2,i] <- lmPCA_mean; sumerrormatrix[2,i] <- sum(Error_lmPCA)
  corLM_mean <- mean(Error_corLM); errormatrix[3,i] <- corLM_mean; sumerrormatrix[3,i] <- sum(Error_corLM)
  stepLM_mean <- mean(Error_stepLM); errormatrix[4,i] <- stepLM_mean; sumerrormatrix[4,i] <- sum(Error_stepLM)
  
  names <- names(errormatrix[1,])
  currentname <- names[i]
  mainname <- paste("Mean Error",currentname)
  mainname2 <- paste("Total Error",currentname)
  
  barplot(errormatrix[,i], names.arg=c("knnPCA","lmPCA","corLM","stepLM"), main = mainname )
  #barplot(sumerrormatrix[,i], names.arg = c("knnPCA","lmPCA","corLM","stepLM"), main = mainname2 )
  # if you want to do the sum of errors.
  
  j <- j + 5
  
}

meanerror <- rowMeans(errormatrix)
barplot(meanerror, main = "Total Mean Error")

par(mfrow=c(1,1))


plot_cumsum(USDmatrix)

USDmatrix <- predictionmatrix[,1:5]

colnames(USDmatrix) <- c("knnPCA","lmPCA","corLM","stepLM","Real Price")
BRLmatrix <- predictionmatrix[,6:10]
colnames(BRLmatrix) <- c("knnPCA","lmPCA","corLM","stepLM","Real Price")
AUDmatrix <- predictionmatrix[,11:15]
colnames(AUDmatrix) <- c("knnPCA","lmPCA","corLM","stepLM","Real Price")
RUBmatrix <- predictionmatrix[,16:20]
colnames(RUBmatrix) <- c("knnPCA","lmPCA","corLM","stepLM","Real Price")
ARSmatrix <- predictionmatrix[,21:25]
colnames(ARSmatrix) <- c("knnPCA","lmPCA","corLM","stepLM","Real Price")
Oilmatrix <- predictionmatrix[,26:30]
colnames(Oilmatrix) <- c("knnPCA","lmPCA","corLM","stepLM","Real Price")


par(mfrow=c(1,1))

plot_cumsum(USDmatrix, "USD",90,120,"topright")
plot_cumsum(BRLmatrix,"BRL",68,110,"topleft")
plot_cumsum(AUDmatrix, "AUD",85,115,"topleft")
plot_cumsum(RUBmatrix,"RUB",70,140,"bottomright")
plot_cumsum(ARSmatrix,"ARS",65,120,"topleft")
plot_cumsum(Oilmatrix,"Oil",55,125,"bottomright")


