#OUTLINE-----------------------------------------------
#Function 1: Data generation for simulation Sudy
#Function 2: Estimators based on Penalty functions
#Function 3: Calculation of evaluation metrics
#------------------------------------------------------
library(glmnet)
library(ncvreg)
#FUNCTION 1: DATA GENERATION---------------------------
simdata <- function(n,CL,k,sc,dist){
  #n : Sample size
  #CL: Censoring level
  #k : No covariates in parametric component should be >=15
  #sc: Selection of scenario, sc=1; \beta=c(5,-3), sc=2; \beta=c(1,-0.5)
  #dist: distribution of covariates dist=1: Normal, dist=2; Uniform

library(Rlab)
library(pracma)
library(psych)
  
  
# Generation of parametric covariates------------------
if (dist==1){
  X <- matrix(rnorm(n*k),nrow=n,ncol=k)  
}
if (dist==2){
  X <- matrix(runif(n*k),nrow=n,ncol=k) 
}
#------------------------------------------------------
#Generation of Regression coefficients-----------------
beta <- 0
if(sc==1){
  for (j in 1:k){
    beta[j] <- 5*I(j>=1 && j<=5)-3*I(j>10 && j<=15)
  }
}
if(sc==2){
  for (j in 1:k){
    beta[j] <- 1*I(j>=1 && j<=5)-0.5*I(j>10 && j<=15)
  }
} 
#Generation of nonparametric covariate and function--------
t <- 0
for (j in 1:n){
  t[j] <- 2.4*((j-0.5)/n)
}
f <- -t*sin(-t^2)
#----------------------------------------------------
#Generation of completely observed PLM---------------
e <- rnorm(n,sd=0.5)
y <- X%*%beta+f+e

#Generation of the right-censored data---------------
delta <-1-rbern(n,CL)     
c <- matrix(0,n,1)
z <- 0
for (i in 1:n){
  if (delta[i]==0){
    while (y[i]<=c[i]){
      c[i]<-rnorm(1,mean(y),sd=sd(y))
    }
  }
  else{
      c[i]<-y[i]
  }
}                          
for (j in 1:n){
  if (y[j]<=c[j]){
    z[j]<-y[j]
  }
  else{
    z[j]<-c[j]
  }
}          
syndata<-function(z,delta){
  library(pracma)
  # where y: right-censored observations, delta: censorship indicator
  n<-length(y)
  M<-0
  yg<-0
  M<-0
  delta2<-0
  #Synthetic data transformation
  y1<-cbind(y,delta)
  y2<-sortrows(y1,k=2)
  delta2<-y2[,2]
  delta2[n]<-1
  sy<-sort(y)
  for (i1 in 1:n){
    Ma<-1
    for (i2 in 1:n){
      mGt<-((n-i2)/(n-i2+1))^((delta2[i2]==0)*(sy[i2]<=sy[i1]))
      Ma<-Ma*mGt
    }
    M[i1]=Ma
    yg[i1]<-(delta[i1]*y[i1])/M[i1]
  }
  return(yg)
}
syn.z <- syndata(z,delta)

out <- new.env()
out$X <- X
out$beta <- beta
out$t <- t
out$f <- f
out$y <- y
out$z <- z
out$c <- c
out$syn.z <- syn.z
out$delta <- delta
out$k <- k
out$CL <- CL
out$n <- n
return(out)
}
#Function which gives descriptive plots for the generated data-----------------
descplots <- function(out){
#out: output argument obtained from the simdata function

#Figure 1: synthetic data, completely observed data and right-censored data----
library(ggplot2)
  library(gridExtra)
t <- out$t
y <- out$y
z <- out$z
syn.z <- out$syn.z
delta <- out$delta
n <- length(y)
group <- ifelse(delta==0, "Right-Censored","Non-Censored")
group2 <- ifelse(delta==0, "Completely Observed","Non-Censored")
group3 <- ifelse(delta==0, "Synthetic (Censored(0))","Synthetic (Observed)")

dfco  <- data.frame(t,y,group2)            #completely observed
dfio  <- data.frame(t,z,group)            #incomplete obs.  
dfsyn <- data.frame(t,syn.z,group3)        #Synthetic responses

p1 <- ggplot()+geom_point(data=dfio, aes( t, z, color = group, shape = group),size = 3)+geom_point(data=dfco, aes(t, y, color = group2, shape = group2),size = 3)+theme(legend.title= element_blank())+ylab("z & y")+ggtitle("(a) Completely (y) and incompletely (z) observed data")
p2 <- ggplot()+geom_point(data=dfio, aes( t, z, color = group, shape = group),size = 3)+geom_point(data=dfsyn, aes( t, syn.z, color = group3, shape = group3),size = 3,stroke=2)+theme(legend.title= element_blank())+ylab("z & yG")+ggtitle("(b) Incomplete (z) and synthetic (yG) observed data")

p3 <- grid.arrange(p1,p2, nrow=1)
#Figure 2: Nonparametric smooth function (f) versus nonp. covariate (t)---------
f <- out$f
fp <- f+rnorm(n)

dff <- data.frame(t,f)
dff2 <- data.frame(t,fp)

p4 <- ggplot()+geom_point(data=dff2,aes(x=t,y=fp),size=2)+geom_line(data=dff,aes(x=t,y=f),size=2)+ylab("f(t)")

outplot <- new.env()

plot(p3)
plot(p4)

outplot$ Censoredplots <- p3
outplot$functionplot <- p4
}

#FUNCTION 2: ESTIMATORS BASED ON DIFFERENT PENALTY FUNCTIONS-------------------
Local_PLM <- function(X,y,t,method,order,nc){
 #method: "Ridge", "Lasso", "SCAD", "MCP", "aLasso", "Enet"
adalasso<-function(X, y,k=10,use.Gram=TRUE,both=TRUE,intercept=TRUE){
    colnames(X)=1:ncol(X)
    n<-length(y)
    cv.adalasso<-NULL
    globalfit<-mylars(X,y,k=k,use.Gram=use.Gram,normalize=TRUE,intercept=intercept)
    coefficients.lasso=globalfit$coefficients
    intercept.lasso=globalfit$intercept
    cv.lasso<-globalfit$cv.lasso
    lambda<-globalfit$lambda
    lambda.lasso<-globalfit$lambda.opt
    coefficients.adalasso=NULL
    lambda.adalasso<-intercept.adalasso<-NULL
    if (use.Gram==TRUE){
      type="covariance"
    }
    if (use.Gram==FALSE){
      type="naive"
    }
    if (both==TRUE){ 
      # cross-validation for adaptive lasso
      all.folds <- split(sample(1:n),rep(1:k,length=n))
      residmat <- matrix(0, length(lambda), k)
      
      for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain<-X[-omit,,drop=FALSE]
        ytrain<-y[-omit]
        Xtest<-X[omit,,drop=FALSE]
        ytest<-y[omit]
        my.lars<-mylars(Xtrain,ytrain,k=k,normalize=TRUE,use.Gram=use.Gram,intercept=intercept)
        coef.lasso<-my.lars$coefficients
        weights <- 1/abs(coef.lasso[ abs(coef.lasso)>0 ])
        #cat(paste("-- non-zero weights ",length(weights),"\n"))
        if (length(weights)==0){
          residmat[,i]<-mean((mean(ytrain)-ytest)^2)
        }
        if (length(weights)==1){
          residmat[,i]=mean((ytest -my.lars$intercept - Xtest%*%coef.lasso)^2)
        }
        if (length(weights)>1){
          XXtrain <- Xtrain[ , names(weights), drop=FALSE]
          XXtest<-Xtest[ , names(weights), drop=FALSE]
          XXtrain <- scale(XXtrain, center=FALSE, scale=weights)
          XXtest<-scale(XXtest, center=FALSE, scale=weights)
          #cat(paste("ncol of XXtrain: ",ncol(XXtrain),"\n"))
          fit<-glmnet(XXtrain,ytrain,type.gaussian=type,standardize=FALSE,intercept=intercept)
          pred<-predict(fit, newx=XXtest, type = "response",s = lambda)
          if (length(omit) == 1){
            pred <- matrix(pred, nrow = 1)
          }
          residmat[, i] <- apply((ytest - pred)^2, 2, mean)
        }
      }
      cv <- apply(residmat, 1, mean)
      cv.adalasso<-min(cv)
      weights <- 1/abs(coefficients.lasso[ abs(coefficients.lasso)>0 ])
      coefficients.adalasso<-rep(0,ncol(X))
      names(coefficients.adalasso)<-1:ncol(X)
      if (length(weights)>0){
        XX <- X[ , names(weights), drop=FALSE]
        if ( length(weights)==1 )  XX <- XX/weights        
        else  XX <- scale(XX, center=FALSE, scale=weights)
        if (length(weights)<=1){
          intercept.adalasso=intercept.lasso 
          coefficients.adalasso<-coefficients.lasso
          lambda.adalasso=0
        }
        else{
          fit<-glmnet(XX,y,type.gaussian=type,standardize=FALSE,intercept=intercept)
          lambda.adalasso<-lambda[which.min(cv)]
          coefficients=predict(fit,type="coefficients",s=lambda.adalasso)
          intercept.adalasso<-coefficients[1]
          coefficients.adalasso[names(weights)]<-coefficients[-1]/weights
        }
      }
    }
    return(list(cv.lasso=cv.lasso,lambda.lasso=lambda.lasso,cv.adalasso=cv.adalasso,lambda.adalasso=lambda.adalasso,intercept.lasso=intercept.lasso, intercept.adalasso=intercept.adalasso, coefficients.lasso=coefficients.lasso,coefficients.adalasso=coefficients.adalasso))
}

mylars<-function (X, y, k = 10,use.Gram=TRUE,normalize=TRUE,intercept=TRUE) 
{
  x<-X
  n<-length(y)
  all.folds <- split(sample(1:n),rep(1:k,length=n))
  
  if (use.Gram==TRUE){
    type="covariance"
  }
  if (use.Gram==FALSE){
    type="naive"
  }
  globalfit<-glmnet(x,y,family="gaussian",standardize=normalize,type.gaussian=type,intercept=intercept)
  lambda<-globalfit$lambda
  residmat <- matrix(0, length(lambda), k)
  for (i in seq(k)) {
    omit <- all.folds[[i]]
    fit <- glmnet(x[-omit, ,drop=FALSE], y[-omit],type.gaussian=type,standardize=normalize,family="gaussian",intercept=intercept)
    fit <- predict(fit, newx=x[omit, , drop = FALSE], type = "response", 
                   s = lambda)
    if (length(omit) == 1) 
      fit <- matrix(fit, nrow = 1)
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
  }
  cv <- apply(residmat, 1, mean)
  cv.lasso<-min(cv)
  cv.error <- sqrt(apply(residmat, 1, var)/k)
  lambda.opt<-lambda[which.min(cv)]
  coefficients=predict(globalfit,type="coefficients",s=lambda.opt)
  inter=coefficients[1]
  coefficients=coefficients[-1]
  names(coefficients)=1:ncol(X)
  object <- list(lambda=lambda,cv=cv,lambda.opt=lambda.opt,cv.lasso=cv.lasso,intercept=inter,coefficients=coefficients)
  invisible(object)
}

aiccfunc<-function(resid,betahat){
    n <- length(y)
    var <- var(resid)
    score<-log(var)+1+((2*(sum(I(betahat!=0))+1))/(n-sum(I(betahat!=0))-2))
    return(score)
  }     #AICc for the tunnig parameter selection
K<-function(W,xx,h){
    t<-(xx-W)/h
    K1<-3/4*(1-t^2)
    return(K1)
  }             #Kernel function for the Local polynomail approach
    dimX <- dim(X)
    k <- dimX[2]
    n<-length(y)
    h <- 0.001
    nolocal <- nc
    L0<-diag(nolocal)
    K1<-matrix(0,n,n)
    WW<-matrix(0,nolocal,nolocal)
    bigS <-matrix(0,n,n)
    bigW <-matrix(0,n,n)
    for (ll in 1:(n-nolocal-1)){
      tt<-seq(t[ll]-1,t[(nolocal+ll-1)]+1,length.out = nolocal)
      Wt <- t[ll:(nolocal+ll-1)]
      yhat<-0
      beta<-0
      for (j in 1:nolocal){
        L0[,j]<-K(Wt[j],tt,h)
        WW[j,j]<-(1/(n*h))*sum(L0)
      }
      ones<-matrix(1,nolocal,1)
      zeros<-matrix(0,nolocal,1)
      ta<-(Wt-tt)
      
      if (order==1){
        xmat<-matrix(c(ones,ta),nolocal,2)
      }
      if (order==2){
        xmat<-matrix(c(ones,ta,ta^2),nolocal,3)
      }
      if (order==3){
        xmat<-matrix(c(ones,ta,ta^2,ta^3),nolocal,4)
      }
      if (order==4){
        xmat<-matrix(c(ones,ta,ta^2,ta^3,ta^4),nolocal,5)
      }
      
      SDL<-xmat%*%solve(t(xmat)%*%WW%*%xmat,tol=1e-100)%*%t(xmat)%*%WW
      bigS[(ll:(nolocal+ll-1)),(ll:(nolocal+ll-1))] <- SDL
      bigW[(ll:(nolocal+ll-1)),(ll:(nolocal+ll-1))] <- WW
    }
    
    
    xtil<-(diag(n)-bigS)%*%X
    ytil<-(diag(n)-bigS)%*%y

  if (method=="Ridge"){
   lam <- seq(0.1,5,length.out = 25)
   AIC.val <- 0
#Choose optimal shrinkage parameter--------------------------------------------
   for (j in 1:25){
     betahatp <- solve(t(xtil)%*%xtil+lam[j]*diag(k))%*%t(xtil)%*%ytil  
     fhatp <- bigS%*%(y-X%*%betahatp)
     yhatp <- X%*%betahatp+fhatp
     resid <- y-yhatp
     AIC.val[j] <- aiccfunc(resid,betahatp)
     if (AIC.val[j]==min(AIC.val)){
       lambda <- lam[j]
   }
   }
   #plot(AIC.val,type="l")
#Obtain PLM ridge-estimator----------------------------------------------------
   betahat <- solve(t(xtil)%*%xtil+lambda*diag(k))%*%t(xtil)%*%ytil  
   fhat <- bigS%*%(y-X%*%betahat)
   yhat <- X%*%betahat+fhat
   resid <- y-yhat
#------------------------------------------------------------------------------
  }
  if (method=="Lasso"){
    lam <- seq(0.01,2,length.out = 25)
    AIC.val <- 0
    for (j in 1:25){
      fit_modelp <- glmnet(xtil, ytil, alpha = 1, lambda = lam[j])
      betahatp <- fit_modelp$beta[,1]
      fhatp <- bigS%*%(y-X%*%betahatp)
      yhatp <- X%*%betahatp+fhatp
      residp <- y-yhatp
      AIC.val[j] <- aiccfunc(residp,betahatp)
      if (AIC.val[j]==min(AIC.val)){
        lambda <- lam[j]
      }
    }
   # plot(AIC.val,type="l")
#Obtain PLM LASSO-estimator----------------------------------------------------
    fit_model <- glmnet(xtil, ytil, alpha = 1, lambda = lambda)
    betahat <- fit_model$beta[,1]
    fhat <- bigS%*%(y-X%*%betahat)
    yhat <- X%*%betahat+fhat
    resid <- y-yhat
  }
  if (method=="aLasso"){
    AIC.val <- 0
    fit_model <- adalasso(xtil,ytil)
    lambda <- fit_model$lambda.adalasso
    betahat <- fit_model$coefficients.adalasso
    fhat <- bigS%*%(y-X%*%betahat)
    yhat <- X%*%betahat+fhat
    resid <- y-yhat
  }
  if (method=="SCAD"){
    lam <- seq(0.001,0.05,length.out = 25)
    AIC.val <- 0
    for (j in 1:25){
      fit_modelp <- ncvreg(xtil, ytil, penalty="SCAD",lambda.min = lam[j])
      betahatp <- as.numeric(fit_modelp$beta[-1,100])
      fhatp <- bigS%*%(y-X%*%betahatp)
      yhatp <- X%*%betahatp+fhatp
      residp <- y-yhatp
      AIC.val[j] <- aiccfunc(residp,betahatp)
      if (AIC.val[j]==min(AIC.val)){
        lambda <- lam[j]
      }
    }
    fit_model <- ncvreg(xtil, ytil, penalty="SCAD",lambda.min = lambda)
    betahat <- as.numeric(fit_model$beta[-1,100])
    fhat <- bigS%*%(y-X%*%betahat)
    yhat <- X%*%betahat+fhat
    resid <- y-yhat
  }
  if (method=="MCP"){
    lam <- seq(0.005,0.01,length.out = 25)
    AIC.val <- 0
    for (j in 1:25){
      fit_modelp <- ncvreg(xtil, ytil, penalty="MCP",lambda.min = lam[j])
      betahatp <- as.numeric(fit_modelp$beta[-1,100])
      fhatp <- bigS%*%(y-X%*%betahatp)
      yhatp <- X%*%betahatp+fhatp
      residp <- y-yhatp
      AIC.val[j] <- aiccfunc(residp,betahatp)
      if (AIC.val[j]==min(AIC.val)){
        lambda <- lam[j]
      }
    }
    fit_model <- ncvreg(xtil, ytil, penalty="MCP",lambda.min = lambda)
    betahat <- as.numeric(fit_model$beta[-1,100])
    fhat <- bigS%*%(y-X%*%betahat)
    yhat <- X%*%betahat+fhat
    resid <- y-yhat
  }
  if (method=="Enet"){
    lam <- seq(0.01,2,length.out = 25)
    alp <- seq(0.01,0.995,length.out = 25)
    AIC.val <- 0
    AIC.val2 <- 0
    for (j in 1:25){
      for (j2 in 1:25){
      fit_modelp <- glmnet(xtil, ytil, alpha = alp[j2], lambda = lam[j])
      betahatp <- fit_modelp$beta[,1]
      fhatp <- bigS%*%(y-X%*%betahatp)
      yhatp <- X%*%betahatp+fhatp
      residp <- y-yhatp
      AIC.val2[j2] <- aiccfunc(residp,betahatp)
      if (AIC.val2[j2]==min(AIC.val2)){
        alpham <- alp[j2]
      }
      }
      AIC.val[j] <- aiccfunc(residp,betahatp)
      if (AIC.val[j]==min(AIC.val)){
        lambda <- lam[j]
      }
    }
    #plot(AIC.val,type="l")
    #Obtain PLM ElasticNet-estimator----------------------------------------------------
    fit_model <- glmnet(xtil, ytil, alpha = alpham, lambda = lambda)
    betahat <- fit_model$beta[,1]
    fhat <- bigS%*%(y-X%*%betahat)
    yhat <- X%*%betahat+fhat
    resid <- y-yhat
  }

out <- new.env()

out$X <- X
out$est.beta <- betahat
out$est.f <- fhat
out$est.y <- yhat
out$chosen.lambda <- lambda
out$AICc <- AIC.val
return(out)
}

#FUNCTION 3: Evaluation Metrics------------------------------------------------
eval.metrics <- function(out,outreal){
#out: obtained object from local_PLM
dimX <- dim(out$X)
y <- outreal$y
k <- dimX[2]
rmse_beta <-sqrt((1/k)*t(out$est.beta-outreal$beta)%*%(out$est.beta-outreal$beta)) #RMSE for BETA
rsq       <- 1-(sum((y-out$est.y)^2)/sum((y-mean(y))^2))                                #R-Square for Model
#Values of Confussion Matrix----------------------------------------------------
a <- 0
b <- 0
c <- 0
d <- 0
for (i in 1:k){
  if (out$est.beta[i]!=0 & outreal$beta[i]!=0){
    d <- d+1
  }
  if (out$est.beta[i]==0 & outreal$beta[i]==0){
    a <- a+1
  }
  if (out$est.beta[i]!=0 & outreal$beta[i]==0){
    c <- c+1
  }
  if (out$est.beta[i]==0 & outreal$beta[i]!=0){
    b <- b+1
  }
}
acc <- (a+d)/(a+b+c+d)
sens <- d/(c+d)
spec <- a/(a+b)
if (spec=="NaN"){
  spec <- 0.01
}
Gscore <- sqrt(sens*spec)

#Evaluation of Nonparametric compoenent-----------------------------------------
mse <- mean((out$est.f-outreal$f)^2)

metrics <- new.env()

metrics$RMSE <- rmse_beta
metrics$Rsq <- rsq
metrics$acc <- acc
metrics$sens <- sens
metrics$spec <- spec
metrics$Gscore <- Gscore
metrics$MSE <- mse

return(metrics)
}

#FUNCTION 4: FIGURES of Simulation Results-------------------------------------
pf_figures <- function(outreal,min,max){
  library(ggplot2)
  adalasso<-function(X, y,k=10,use.Gram=TRUE,both=TRUE,intercept=TRUE,lambdaout){
    colnames(X)=1:ncol(X)
    n<-length(y)
    cv.adalasso<-NULL
    globalfit<-mylars(X,y,k=k,use.Gram=use.Gram,normalize=TRUE,intercept=intercept)
    coefficients.lasso=globalfit$coefficients
    intercept.lasso=globalfit$intercept
    cv.lasso<-globalfit$cv.lasso
    lambda<-globalfit$lambda
    lambda.lasso<-globalfit$lambda.opt
    coefficients.adalasso=NULL
    lambda.adalasso<-intercept.adalasso<-NULL
    if (use.Gram==TRUE){
      type="covariance"
    }
    if (use.Gram==FALSE){
      type="naive"
    }
    if (both==TRUE){ 
      # cross-validation for adaptive lasso
      all.folds <- split(sample(1:n),rep(1:k,length=n))
      residmat <- matrix(0, length(lambda), k)
      
      for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain<-X[-omit,,drop=FALSE]
        ytrain<-y[-omit]
        Xtest<-X[omit,,drop=FALSE]
        ytest<-y[omit]
        my.lars<-mylars(Xtrain,ytrain,k=k,normalize=TRUE,use.Gram=use.Gram,intercept=intercept)
        coef.lasso<-my.lars$coefficients
        weights <- 1/abs(coef.lasso[ abs(coef.lasso)>0 ])
        #cat(paste("-- non-zero weights ",length(weights),"\n"))
        if (length(weights)==0){
          residmat[,i]<-mean((mean(ytrain)-ytest)^2)
        }
        if (length(weights)==1){
          residmat[,i]=mean((ytest -my.lars$intercept - Xtest%*%coef.lasso)^2)
        }
        if (length(weights)>1){
          XXtrain <- Xtrain[ , names(weights), drop=FALSE]
          XXtest<-Xtest[ , names(weights), drop=FALSE]
          XXtrain <- scale(XXtrain, center=FALSE, scale=weights)
          XXtest<-scale(XXtest, center=FALSE, scale=weights)
          #cat(paste("ncol of XXtrain: ",ncol(XXtrain),"\n"))
          fit<-glmnet(XXtrain,ytrain,type.gaussian=type,standardize=FALSE,intercept=intercept)
          pred<-predict(fit, newx=XXtest, type = "response",s = lambda)
          if (length(omit) == 1){
            pred <- matrix(pred, nrow = 1)
          }
          residmat[, i] <- apply((ytest - pred)^2, 2, mean)
        }
      }
      cv <- apply(residmat, 1, mean)
      cv.adalasso<-min(cv)
      weights <- 1/abs(coefficients.lasso[ abs(coefficients.lasso)>0 ])
      coefficients.adalasso<-rep(0,ncol(X))
      names(coefficients.adalasso)<-1:ncol(X)
      if (length(weights)>0){
        XX <- X[ , names(weights), drop=FALSE]
        if ( length(weights)==1 )  XX <- XX/weights        
        else  XX <- scale(XX, center=FALSE, scale=weights)
        if (length(weights)<=1){
          intercept.adalasso=intercept.lasso 
          coefficients.adalasso<-coefficients.lasso
          lambda.adalasso=0
        }
        else{
          fit<-glmnet(XX,y,type.gaussian=type,standardize=FALSE,intercept=intercept)
          lambda.adalasso<-lambda[which.min(cv)]
          coefficients=predict(fit,type="coefficients",s=lambdaout)
          intercept.adalasso<-coefficients[1]
          coefficients.adalasso[names(weights)]<-coefficients[-1]/weights
        }
      }
    }
    return(list(cv.lasso=cv.lasso,lambda.lasso=lambda.lasso,cv.adalasso=cv.adalasso,lambda.adalasso=lambda.adalasso,intercept.lasso=intercept.lasso, intercept.adalasso=intercept.adalasso, coefficients.lasso=coefficients.lasso,coefficients.adalasso=coefficients.adalasso))
  }
  
  mylars<-function (X, y, k = 10,use.Gram=TRUE,normalize=TRUE,intercept=TRUE) 
  {
    x<-X
    n<-length(y)
    all.folds <- split(sample(1:n),rep(1:k,length=n))
    
    if (use.Gram==TRUE){
      type="covariance"
    }
    if (use.Gram==FALSE){
      type="naive"
    }
    globalfit<-glmnet(x,y,family="gaussian",standardize=normalize,type.gaussian=type,intercept=intercept)
    lambda<-globalfit$lambda
    residmat <- matrix(0, length(lambda), k)
    for (i in seq(k)) {
      omit <- all.folds[[i]]
      fit <- glmnet(x[-omit, ,drop=FALSE], y[-omit],type.gaussian=type,standardize=normalize,family="gaussian",intercept=intercept)
      fit <- predict(fit, newx=x[omit, , drop = FALSE], type = "response", 
                     s = lambda)
      if (length(omit) == 1) 
        fit <- matrix(fit, nrow = 1)
      residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    }
    cv <- apply(residmat, 1, mean)
    cv.lasso<-min(cv)
    cv.error <- sqrt(apply(residmat, 1, var)/k)
    lambda.opt<-lambda[which.min(cv)]
    coefficients=predict(globalfit,type="coefficients",s=lambda.opt)
    inter=coefficients[1]
    coefficients=coefficients[-1]
    names(coefficients)=1:ncol(X)
    object <- list(lambda=lambda,cv=cv,lambda.opt=lambda.opt,cv.lasso=cv.lasso,intercept=inter,coefficients=coefficients)
    invisible(object)
  }
  X <- outreal$X
  y <- outreal$syn.z
  t <- outreal$t
  K<-function(W,xx,h){
    t<-(xx-W)/h
    K1<-3/4*(1-t^2)
    return(K1)
  }             #Kernel function for the Local polynomail approach
  h <- 0.1
  order <- 2
  dimX <- dim(X)
  k <- dimX[2]
  n<-length(y)
  nolocal <- 7
  L0<-diag(nolocal)
  K1<-matrix(0,n,n)
  WW<-matrix(0,nolocal,nolocal)
  bigS <-matrix(0,n,n)
  bigW <-matrix(0,n,n)
  for (ll in 1:(n-6)){
    tt<-seq(t[ll]-1,t[(nolocal+ll-1)]+1,length.out = nolocal)
    Wt <- t[ll:(nolocal+ll-1)]
    yhat<-0
    beta<-0
    for (j in 1:nolocal){
      L0[,j]<-K(Wt[j],tt,h)
      WW[j,j]<-(1/(n*h))*sum(L0)
    }
    ones<-matrix(1,nolocal,1)
    zeros<-matrix(0,nolocal,1)
    ta<-(Wt-tt)
    
    if (order==1){
      xmat<-matrix(c(ones,ta),nolocal,2)
    }
    if (order==2){
      xmat<-matrix(c(ones,ta,ta^2),nolocal,3)
    }
    if (order==3){
      xmat<-matrix(c(ones,ta,ta^2,ta^3),nolocal,4)
    }
    if (order==4){
      xmat<-matrix(c(ones,ta,ta^2,ta^3,ta^4),nolocal,5)
    }
    
    SDL<-xmat%*%solve(t(xmat)%*%WW%*%xmat,tol=1e-100)%*%t(xmat)%*%WW
    bigS[(ll:(nolocal+ll-1)),(ll:(nolocal+ll-1))] <- SDL
    bigW[(ll:(nolocal+ll-1)),(ll:(nolocal+ll-1))] <- WW
  }
  xtil<-(diag(n)-bigS)%*%X
  ytil<-(diag(n)-bigS)%*%y
  
  betahat_ridge <- matrix(0,k,50)
  betahat_ridgex <- matrix(0,k,50)
  betahat_lasso <- matrix(0,k,50)
  betahat_alasso<- matrix(0,k,50)
  betahat_SCAD  <- matrix(0,k,50)
  betahat_MCP   <- matrix(0,k,50)
  betahat_Enet  <- matrix(0,k,50)
  #RIDGE
  lambdaseq <- seq(min,max,length.out = 50)
  for (i in 1:50){
  betahat_ridge[,i] <- solve(t(xtil)%*%xtil+lambdaseq[i]*diag(k))%*%t(xtil)%*%ytil 
  #LASSO
  fit_model <- glmnet(xtil, ytil, alpha = 1, lambda = lambdaseq[i])
  betahat_lasso[,i] <- fit_model$beta[,1]
  #ADAPTIVE-LASSO
  fit_model <- adalasso(xtil,ytil,lambdaout =lambdaseq[i])
  lambda <- fit_model$lambda.adalasso
  betahat_alasso[,i] <- fit_model$coefficients.adalasso
  #SCAD
  fit_model        <- ncvreg(xtil, ytil, penalty="SCAD",lambda.min =lambdaseq[i]/2)
  betahat_SCAD[,i] <- as.numeric(fit_model$beta[-1,100])
  #MCP
  fit_model        <- ncvreg(xtil, ytil, penalty="MCP",lambda.min =lambdaseq[i]/2)
  betahat_MCP[,i]  <- as.numeric(fit_model$beta[-1,100])
  #ElasticNET
  fit_model        <- glmnet(xtil, ytil, alpha = 0.7, lambda = lambdaseq[i])
  betahat_Enet[,i] <- fit_model$beta[,1]
  }
  library(reshape)
  df_ridge <- data.frame(lambdaseq,t(betahat_ridge))
  df_ridge <- melt(df_ridge, id.vars = "lambdaseq")
  
  df_lasso <- data.frame(lambdaseq,t(betahat_lasso))
  df_lasso <- melt(df_lasso, id.vars = "lambdaseq")
  
  df_alasso <- data.frame(lambdaseq,t(betahat_alasso))
  df_alasso <- melt(df_alasso, id.vars = "lambdaseq")
  
  df_scad <- data.frame(lambdaseq,t(betahat_SCAD))
  df_scad <- melt(df_scad, id.vars = "lambdaseq")
  
  df_mcp <- data.frame(lambdaseq,t(betahat_MCP))
  df_mcp <- melt(df_mcp, id.vars = "lambdaseq")
  
  df_enet <- data.frame(lambdaseq,t(betahat_Enet))
  df_enet <- melt(df_enet, id.vars = "lambdaseq")
  #-------------------------------------------------------------------------------
  pridge <- ggplot(df_ridge, aes(x = lambdaseq, y = value, colour = variable)) + geom_line(lwd=1)+theme(legend.position="none")+ylab("Beta")+xlab("lambda")+ggtitle("(a) Ridge")+geom_hline(yintercept=0)
  plasso <- ggplot(df_lasso, aes(x = lambdaseq, y = value, colour = variable)) + geom_line(lwd=1)+theme(legend.position="none")+ylab("Beta")+xlab("lambda")+ggtitle("(b) Lasso")+geom_hline(yintercept=0)
  palasso <- ggplot(df_alasso, aes(x = lambdaseq, y = value, colour = variable)) + geom_line(lwd=1)+theme(legend.position="none")+ylab("Beta")+xlab("lambda")+ggtitle("(c) aLasso")+geom_hline(yintercept=0)
  pscad <- ggplot(df_scad, aes(x = lambdaseq, y = value, colour = variable)) + geom_line(lwd=1)+theme(legend.position="none")+ylab("Beta")+xlab("lambda")+ggtitle("(d) SCAD")+geom_hline(yintercept=0)
  pmcp <- ggplot(df_mcp, aes(x = lambdaseq, y = value, colour = variable)) + geom_line(lwd=1)+theme(legend.position="none")+ylab("Beta")+xlab("lambda")+ggtitle("(e) MCP")+geom_hline(yintercept=0)
  penet <- ggplot(df_enet, aes(x = lambdaseq, y = value, colour = variable)) + geom_line(lwd=1)+theme(legend.position="none")+ylab("Beta")+xlab("lambda")+ggtitle("(f) ElasticNet")+geom_hline(yintercept=0)
  
  
  library(gridExtra)
  grid.arrange(pridge,plasso,palasso,pscad,pmcp,penet,nrow=1)
  
}


#sim.fig <- function(out,outreal){
  
#}
  

