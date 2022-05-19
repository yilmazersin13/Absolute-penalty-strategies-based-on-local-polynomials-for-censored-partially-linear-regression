#SIMUALTION STUDY FOR lOCAL-PLM WITH PENALTY FUNCTIONS-------------------------
#s1: Loop index for sample size, n=50,150,300
#s2: Loop index for CL=10%, 30%
#s3: Loop index for k=15,40
sim <- 2
a <- 1 #Combination index (upto 24)
#EMPTIES-----------------------------------------------------------------------
RMSE <- matrix(0,sim,10)
colnames(RMSE) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
RMSEx <- list(matrix(0,sim,10))
MSE <- matrix(0,sim,10)
colnames(MSE) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
MSEx <- list(matrix(0,sim,10))
Rsquare <- matrix(0,sim,10)
colnames(Rsquare) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
Rsquarex <- list(matrix(0,sim,10))
Gscore <- matrix(0,sim,10)
colnames(Gscore) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
Gscorex <- list(matrix(0,sim,10))
Sens <- matrix(0,sim,10)
colnames(Sens) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
Sensx <- list(matrix(0,sim,10))
Spec <- matrix(0,sim,10)
colnames(Spec) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
Specx <- list(matrix(0,sim,10))
Accur <- matrix(0,sim,10)
colnames(Accur) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scenario","Sample Size","CL","No. Covariate")
Accurx <- list(matrix(0,sim,10))

Beta_Ridgex <- list()
Beta_Lassox <- list()
Beta_aLassox <- list()
Beta_SCADx <- list()
Beta_MCPx <- list()
Beta_Enetx <- list()
f_Ridgex <- list()
f_Lassox <- list()
f_aLassox <- list()
f_SCADx <- list()
f_MCPx <- list()
f_Enetx <- list()
y_Ridgex <- list()
y_Lassox <- list()
y_aLassox <- list()
y_SCADx <- list()
y_MCPx <- list()
y_Enetx <- list()
#-------------------------------------------------------------------------------
for (s1 in 1:3){
  for (s2 in 1:2){
    for (s3 in 1:2){
      for(sc in 1:2){
      if (s1==1){n=50} 
      if (s1==2){n=150} 
      if (s1==3){n=300}
      CL <- ifelse(s2==1,0.10,0.30)
      k <- ifelse(s3==1,15,40)
#EMPTY BETA & f(z)-------------------------------------------------------------
      Est.Beta.Ridge  <- matrix(0,k,sim)
      Est.Beta.Lasso  <- matrix(0,k,sim)
      Est.Beta.aLasso <- matrix(0,k,sim)
      Est.Beta.SCAD   <- matrix(0,k,sim)
      Est.Beta.MCP    <- matrix(0,k,sim)
      Est.Beta.Enet   <- matrix(0,k,sim)
      
      Est.f.Ridge  <- matrix(0,n,sim)
      Est.f.Lasso  <- matrix(0,n,sim)
      Est.f.aLasso <- matrix(0,n,sim)
      Est.f.SCAD   <- matrix(0,n,sim)
      Est.f.MCP    <- matrix(0,n,sim)
      Est.f.Enet   <- matrix(0,n,sim)
      
      Est.y.Ridge  <- matrix(0,n,sim)
      Est.y.Lasso  <- matrix(0,n,sim)
      Est.y.aLasso <- matrix(0,n,sim)
      Est.y.SCAD   <- matrix(0,n,sim)
      Est.y.MCP    <- matrix(0,n,sim)
      Est.y.Enet   <- matrix(0,n,sim)
#------------------------------------------------------------------------------
        for (s in 1:sim){
      data <- simdata(n,CL,k,sc,1)
      #descplots(data)
      X <- data$X
      yg <- data$syn.z
      t <- data$t
      order <- 2
      #ESTIMATED OBJECTS BASED ON SIX METHODS-----------------------------------
      obj_Ridge <- Local_PLM(X,yg,t,"Ridge",order,n/5)
      obj_Lasso <- Local_PLM(X,yg,t,"Lasso",order,n/5)
      obj_aLasso <- Local_PLM(X,yg,t,"aLasso",order,n/5)
      obj_SCAD <- Local_PLM(X,yg,t,"SCAD",order,n/5)
      obj_MCP <- Local_PLM(X,yg,t,"MCP",order,n/5)
      obj_Enet <- Local_PLM(X,yg,t,"Enet",order,n/5)
      
      Metrics_Ridge <- eval.metrics(obj_Ridge,data)
      Metrics_Lasso <- eval.metrics(obj_Lasso,data)
      Metrics_aLasso <- eval.metrics(obj_aLasso,data)
      Metrics_SCAD <- eval.metrics(obj_SCAD,data)
      Metrics_MCP <- eval.metrics(obj_MCP,data)
      Metrics_Enet <- eval.metrics(obj_Enet,data)
      
     
      RMSE[s,] <- c(Metrics_Ridge$RMSE,Metrics_Lasso$RMSE,Metrics_aLasso$RMSE,Metrics_SCAD$RMSE,Metrics_MCP$RMSE,Metrics_Enet$RMSE,sc,n,CL,k) 
      Rsquare[s,] <- c(Metrics_Ridge$Rsq,Metrics_Lasso$Rsq,Metrics_aLasso$Rsq,Metrics_SCAD$Rsq,Metrics_MCP$Rsq,Metrics_Enet$Rsq,sc,n,CL,k) 
      MSE[s,] <- c(Metrics_Ridge$MSE,Metrics_Lasso$MSE,Metrics_aLasso$MSE,Metrics_SCAD$MSE,Metrics_MCP$MSE,Metrics_Enet$MSE,sc,n,CL,k) 
      Gscore[s,] <- c(Metrics_Ridge$Gscore,Metrics_Lasso$Gscore,Metrics_aLasso$Gscore,Metrics_SCAD$Gscore,Metrics_MCP$Gscore,Metrics_Enet$Gscore,sc,n,CL,k) 
      Sens[s,] <-  c(Metrics_Ridge$sens,Metrics_Lasso$sens,Metrics_aLasso$sens,Metrics_SCAD$sens,Metrics_MCP$sens,Metrics_Enet$sens,sc,n,CL,k) 
      Spec[s,] <-  c(Metrics_Ridge$spec,Metrics_Lasso$spec,Metrics_aLasso$spec,Metrics_SCAD$spec,Metrics_MCP$spec,Metrics_Enet$spec,sc,n,CL,k) 
      Accur[s,] <-  c(Metrics_Ridge$acc,Metrics_Lasso$acc,Metrics_aLasso$acc,Metrics_SCAD$acc,Metrics_MCP$acc,Metrics_Enet$acc,sc,n,CL,k) 
      
      Est.Beta.Ridge[,s]  <- obj_Ridge$est.beta
      Est.Beta.Lasso[,s]  <- obj_Lasso$est.beta
      Est.Beta.aLasso[,s] <- obj_aLasso$est.beta
      Est.Beta.SCAD[,s]   <- obj_SCAD$est.beta
      Est.Beta.MCP[,s]    <- obj_MCP$est.beta
      Est.Beta.Enet[,s]   <- obj_Enet$est.beta
      
      Est.f.Ridge[,s]  <- obj_Ridge$est.f
      Est.f.Lasso[,s]  <- obj_Lasso$est.f
      Est.f.aLasso[,s] <- obj_aLasso$est.f
      Est.f.SCAD[,s]   <- obj_SCAD$est.f
      Est.f.MCP[,s]    <- obj_MCP$est.f
      Est.f.Enet[,s]   <- obj_Enet$est.f
      
      Est.y.Ridge[,s]  <- obj_Ridge$est.y
      Est.y.Lasso[,s]  <- obj_Lasso$est.y
      Est.y.aLasso[,s] <- obj_aLasso$est.y
      Est.y.SCAD[,s]   <- obj_SCAD$est.y
      Est.y.MCP[,s]    <- obj_MCP$est.y
      Est.y.Enet[,s]   <- obj_Enet$est.y
      
      message("Simulation ends",s)
      }
      RMSEx[[a]]<-RMSE
      Rsquarex[[a]]<-Rsquare
      MSEx[[a]]<-MSE
      Gscorex[[a]]<-Gscore
      Sensx[[a]]<-Sens
      Specx[[a]]<-Spec
      Accurx[[a]]<-Accur
      
      Beta_Ridgex[[a]] <- Est.Beta.Ridge
      Beta_Lassox[[a]] <- Est.Beta.Lasso
      Beta_aLassox[[a]] <- Est.Beta.aLasso
      Beta_SCADx[[a]] <- Est.Beta.SCAD
      Beta_MCPx[[a]] <- Est.Beta.MCP
      Beta_Enetx[[a]] <- Est.Beta.Enet
      
      f_Ridgex[[a]] <- Est.f.Ridge
      f_Lassox[[a]] <- Est.f.Lasso
      f_aLassox[[a]] <- Est.f.aLasso
      f_SCADx[[a]] <- Est.f.SCAD
      f_MCPx[[a]] <- Est.f.MCP
      f_Enetx[[a]] <- Est.f.Enet
      
      y_Ridgex[[a]] <- Est.y.Ridge
      y_Lassox[[a]] <- Est.y.Lasso
      y_aLassox[[a]] <- Est.y.aLasso
      y_SCADx[[a]] <- Est.y.SCAD
      y_MCPx[[a]] <- Est.y.MCP
      y_Enetx[[a]] <- Est.y.Enet
    
      a <- a+1
      message("Scenario",sc)
      message("Results for n=",n)
      message("Censoring level=",CL)
      message("no. covariate=",k)
      }
    }
  }
}

