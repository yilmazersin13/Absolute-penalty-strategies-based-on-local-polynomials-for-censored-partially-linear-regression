#DESIGNING SIMULATION RESULTS--------------------------------------------------
#data <- simdata(150,0.30,15,1,1)
#Scenario1, n=50, k=15, CL=0.10
Est.Beta.Lasso  <- Beta_Lassox[[1]]
Est.Beta.aLasso <- Beta_aLassox[[1]]
Est.Beta.SCAD   <- Beta_SCADx[[1]]
Est.Beta.MCP    <- Beta_MCPx[[1]]
Est.Beta.Enet   <- Beta_Enetx[[1]]



#Scenario1, n=50, k=15, CL=0.30
Est.Beta.Lasso  <- Beta_Lassox[[5]]
Est.Beta.aLasso <- Beta_aLassox[[5]]
Est.Beta.SCAD   <- Beta_SCADx[[5]]
Est.Beta.MCP    <- Beta_MCPx[[5]]
Est.Beta.Enet   <- Beta_Enetx[[5]]

#Scenario2, n=150, k=40, CL=0.10
Est.Beta.Lasso  <- Beta_Lassox[[12]]
Est.Beta.aLasso <- Beta_aLassox[[12]]
Est.Beta.SCAD   <- Beta_SCADx[[12]]
Est.Beta.MCP    <- Beta_MCPx[[12]]
Est.Beta.Enet   <- Beta_Enetx[[12]]

#Scenario2, n=300, k=40, CL=0.30
Est.Beta.Lasso  <- Beta_Lassox[[24]]
Est.Beta.aLasso <- Beta_aLassox[[24]]
Est.Beta.SCAD   <- Beta_SCADx[[24]]
Est.Beta.MCP    <- Beta_MCPx[[24]]
Est.Beta.Enet   <- Beta_Enetx[[24]]

#CONSTRUCTION OF VARIABLES FOR BARPLOTS------------------------------------------
b.lasso <- 0
b.alasso <- 0
b.scad <- 0
b.mcp <- 0
b.enet <- 0
k <- 40
  for (j in 1:k){
    nz.lasso <- 0
    z.lasso <- 0 
    nz.alasso <- 0
    z.alasso <- 0 
    nz.scad <- 0
    z.scad <- 0 
    nz.mcp <- 0
    z.mcp <- 0 
    nz.enet <- 0
    z.enet <- 0 
    
    for (i in 1:sim){
      if (data$beta[j]!=0 & Est.Beta.Lasso[j,i]!=0){
        nz.lasso <- nz.lasso+1
        b.lasso[j] <- nz.lasso-1
      }  
      if (data$beta[j]==0 & Est.Beta.Lasso[j,i]==0){
        z.lasso <- z.lasso+1
        b.lasso[j] <- z.lasso-1
      }  
      if (data$beta[j]!=0 & Est.Beta.aLasso[j,i]!=0){
        nz.alasso <- nz.alasso+1
        b.alasso[j] <- nz.alasso-1
      }  
      if (data$beta[j]==0 & Est.Beta.aLasso[j,i]==0){
        z.alasso <- z.alasso+1
        b.alasso[j] <- z.alasso-1
      }  
      if (data$beta[j]!=0 & Est.Beta.SCAD[j,i]!=0){
        nz.scad <- nz.scad+1
        b.scad[j] <- nz.scad-1
      }  
      if (data$beta[j]==0 & Est.Beta.SCAD[j,i]==0){
        z.scad <- z.scad+1
        b.scad[j] <- z.scad-1
      }  
      if (data$beta[j]!=0 & Est.Beta.MCP[j,i]!=0){
        nz.mcp <- nz.mcp+1
        b.mcp[j] <- nz.mcp-1
      }  
      if (data$beta[j]==0 & Est.Beta.MCP[j,i]==0){
        z.mcp <- z.mcp+1
        b.mcp[j] <- z.mcp-1
      } 
      if (data$beta[j]!=0 & Est.Beta.Enet[j,i]!=0){
        nz.enet <- nz.enet+1
        b.enet[j] <- nz.enet-1
      }  
      if (data$beta[j]==0 & Est.Beta.Enet[j,i]==0){
        z.enet <- z.enet+1
        b.enet[j] <- z.enet-1
      }  
    }
  }
 b.enet[40] <- NA
 #Selection frequencies of true betas (%)---------------------------------------
 library(dplyr)
 library(ggplot2)
 Results <- data.frame(b.lasso/sim,b.alasso/sim,b.scad/sim,b.mcp/sim,b.enet/sim)

p.lasso <- b.lasso/sim
p.alasso <- b.alasso/sim
p.scad <- b.scad/sim
p.mcp <- b.mcp/sim
p.enet <- b.enet/sim
index <- c(1:k)

dflasso <- data.frame(index,p.lasso)
dfalasso <- data.frame(index,p.alasso)
dfscad <- data.frame(index,p.scad)
dfmcp <- data.frame(index,p.mcp)
dfenet <- data.frame(index,p.enet)
 
dflasso <- dflasso %>%
  mutate(cond = case_when(
    index<=5  ~ 'black',
    index>10 & index<=15 ~ 'black',
    TRUE ~ 'gray'   #anything that does not meet the criteria above
  ))
dfalasso <- dfalasso %>%
  mutate(cond = case_when(
    index<=5  ~ 'black',
    index>10 & index<=15 ~ 'black',
    TRUE ~ 'gray'   #anything that does not meet the criteria above
  ))
dfscad <- dfscad %>%
  mutate(cond = case_when(
    index<=5  ~ 'black',
    index>10 & index<=15 ~ 'black',
    TRUE ~ 'gray'   #anything that does not meet the criteria above
  ))
dfmcp <- dfmcp %>%
  mutate(cond = case_when(
    index<=5  ~ 'black',
    index>10 & index<=15 ~ 'black',
    TRUE ~ 'gray'   #anything that does not meet the criteria above
  ))
dfenet <- dfenet %>%
  mutate(cond = case_when(
    index<=5  ~ 'black',
    index>10 & index<=15 ~ 'black',
    TRUE ~ 'gray'   #anything that does not meet the criteria above
  ))

pll <- ggplot(data = dflasso, aes(x = index, y =p.lasso)) +
  geom_bar(stat = "identity", color = 'black', aes(fill = cond)) +
  scale_fill_identity()+xlab("No.Coefficient")+ylab("Selection frequency (%)")+ggtitle("(a) Modified Lasso (LL)")
pall <- ggplot(data = dfalasso, aes(x = index, y =p.alasso)) +
  geom_bar(stat = "identity", color = 'black', aes(fill = cond)) +
  scale_fill_identity()+xlab("No.Coefficient")+ylab("Selection frequency (%)")+ggtitle("(b) Modified aLasso (aLL)")
psl <- ggplot(data = dfscad, aes(x = index, y =p.scad)) +
  geom_bar(stat = "identity", color = 'black', aes(fill = cond)) +
  scale_fill_identity()+xlab("No.Coefficient")+ylab("Selection frequency (%)")+ggtitle("(c) Modified SCAD (SL)")
pmcl <- ggplot(data = dfmcp, aes(x = index, y =p.mcp)) +
  geom_bar(stat = "identity", color = 'black', aes(fill = cond)) +
  scale_fill_identity()+xlab("No.Coefficient")+ylab("Selection frequency (%)")+ggtitle("(d) Modified MCP (MCL)")
penet <- ggplot(data = dfenet, aes(x = index, y =p.enet)) +
  geom_bar(stat = "identity", color = 'black', aes(fill = cond)) +
  scale_fill_identity()+xlab("No.Coefficient")+ylab("Selection frequency (%)")+ggtitle("(e) Modified Elasticnet (ENL)")

gridExtra::grid.arrange(pll,pall,psl,pmcl,penet, 
                        layout_matrix=rbind(c(1,1,2,2),c(3,3,4,4),c(NA,5,5,NA)))
#------------------------------------------------------------------------------
Gtable <- matrix(0,24,10)
Acctable <- matrix(0,24,10)
RMSEtable <- matrix(0,24,10)
Rsqtable <- matrix(0,24,10)
#PLOTS and TABLES FOR EVAL METRICS--------------------------------------------
for (i in 1:24){
Gtable[i,]    <- colMeans(Gscorex[[i]]) 
Acctable[i,]  <- colMeans(Accurx[[i]])
RMSEtable[i,] <- colMeans(RMSEx[[i]])
Rsqtable[i,]  <- colMeans(Rsquarex[[i]])
}
write.table(Gtable,"Gscores.txt")
write.table(Acctable,"Accuracy.txt")
write.table(RMSEtable,"RMSEbeta.txt")
write.table(Rsqtable,"Rsquare.txt")

#RESULTS of NONPARAMETRIC REGRESSION-------------------------------------------
#Scenario 1, k=40, CL=10%, n=150
n <- 150
f_Ridge <- f_Ridgex[[11]]
f_Lasso <- f_Lassox[[11]]
f_aLasso <- f_aLassox[[11]]
f_SCAD <- f_SCADx[[11]]
f_MCP <- f_MCPx[[11]]
f_Enet <- f_Enetx[[11]]
t <- 0
for (i in 1:n){t[i] <- 2.4*(i-0.5)/n}
f <- -t*sin(-t^2)
#Scenario 1, k=40, CL=30%, n=150
n <- 150
f_Ridge <- f_Ridgex[[15]]
f_Lasso <- f_Lassox[[15]]
f_aLasso <- f_aLassox[[15]]
f_SCAD <- f_SCADx[[15]]
f_MCP <- f_MCPx[[15]]
f_Enet <- f_Enetx[[15]]
t <- 0
for (i in 1:n){t[i] <- 2.4*(i-0.5)/n}
f <- -t*sin(-t^2)
#Scenario 2, k=15, CL=10%, n=50
n <- 300
f_Ridge <- f_Ridgex[[18]]
f_Lasso <- f_Lassox[[18]]
f_aLasso <- f_aLassox[[18]]
f_SCAD <- f_SCADx[[18]]
f_MCP <- f_MCPx[[18]]
f_Enet <- f_Enetx[[18]]
t <- 0
for (i in 1:n){t[i] <- 2.4*(i-0.5)/n}
f <- -t*sin(-t^2)
#Scenario 2, k=15, CL=30%, n=300
n <- 300
f_Ridge <- f_Ridgex[[22]]
f_Lasso <- f_Lassox[[22]]
f_aLasso <- f_aLassox[[22]]
f_SCAD <- f_SCADx[[22]]
f_MCP <- f_MCPx[[22]]
f_Enet <- f_Enetx[[22]]
t <- 0
for (i in 1:n){t[i] <- 2.4*(i-0.5)/n}
f <- -t*sin(-t^2)
#------------------------------------------------------------------------------
df_real <- data.frame(t,f)
df_lasso <- data.frame(t,(f_Lasso))
df_lasso <- melt(df_lasso, id.vars = "t")

df_ridge <- data.frame(t,(f_Ridge))
df_ridge <- melt(df_ridge, id.vars = "t")

df_alasso <- data.frame(t,(f_aLasso))
df_alasso <- melt(df_alasso, id.vars = "t")

df_scad <- data.frame(t,(f_SCAD))
df_scad <- melt(df_scad, id.vars = "t")

df_mcp <- data.frame(t,(f_MCP))
df_mcp <- melt(df_mcp, id.vars = "t")

df_enet <- data.frame(t,(f_Enet))
df_enet <- melt(df_enet, id.vars = "t")
#-------------------------------------------------------------------------------
plasso <- ggplot() + geom_line(data=df_lasso, aes(x = t, y = value/2))+geom_line(data=df_real, aes(x = t, y = f, colour = "black"),lwd=1)+theme(legend.position="none")+ylab("f(lasso)")+xlab("t")+ggtitle("(a) LL estimator")
pridge <- ggplot() + geom_line(data=df_ridge, aes(x = t, y = value/2))+geom_line(data=df_real, aes(x = t, y = f, colour = "black"),lwd=1)+theme(legend.position="none")+ylab("f(Ridge)")+xlab("t")+ggtitle("(b) RL estimator")
palasso <- ggplot() + geom_line(data=df_alasso, aes(x = t, y = value/2))+geom_line(data=df_real, aes(x = t, y = f, colour = "black"),lwd=1)+theme(legend.position="none")+ylab("f(aLasso)")+xlab("t")+ggtitle("(c) aLL estimator")
pscad <- ggplot() + geom_line(data=df_scad, aes(x = t, y = value/2))+geom_line(data=df_real, aes(x = t, y = f, colour = "black"),lwd=1)+theme(legend.position="none")+ylab("f(Scad)")+xlab("t")+ggtitle("(d) SL estimator")
pmcp <- ggplot() + geom_line(data=df_mcp, aes(x = t, y = value/2))+geom_line(data=df_real, aes(x = t, y = f, colour = "black"),lwd=1)+theme(legend.position="none")+ylab("f(MCP)")+xlab("t")+ggtitle("(e) MCL estimator")
penet <- ggplot() + geom_line(data=df_enet, aes(x = t, y = value/2))+geom_line(data=df_real, aes(x = t, y = f, colour = "black"),lwd=1)+theme(legend.position="none")+ylab("f(Enet)")+xlab("t")+ggtitle("(f) ENL estimator")

grid.arrange(plasso,pridge,palasso,pscad,pmcp,penet,nrow=3)

#INDIVIDUAL PLOTS FOR NP FUNCTIONS----------------------------------------------------------------------------------------------
#RESULTS of NONPARAMETRIC REGRESSION-------------------------------------------
#Scenario 1, k=15, CL=10%, n=50
n <- 50
f_Ridge1 <- rowMeans(f_Ridgex[[1]])  
f_Lasso1 <- rowMeans(f_Lassox[[1]])
f_aLasso1 <- rowMeans(f_aLassox[[1]])
f_SCAD1 <- rowMeans(f_SCADx[[1]])
f_MCP1 <- rowMeans(f_MCPx[[1]])
f_Enet1 <- rowMeans(f_Enetx[[1]])
t1 <- 0
for (i in 1:n){t1[i] <- 2.4*(i-0.5)/n}
f1 <- -t1*sin(-t1^2)
#Scenario 1, k=40, CL=30%, n=150
n <- 150
f_Ridge2 <- rowMeans(f_Ridgex[[15]])
f_Lasso2 <- rowMeans(f_Lassox[[15]])
f_aLasso2 <- rowMeans(f_aLassox[[15]])
f_SCAD2 <- rowMeans(f_SCADx[[15]])
f_MCP2 <- rowMeans(f_MCPx[[15]])
f_Enet2 <- rowMeans(f_Enetx[[15]])
t2 <- 0
for (i in 1:n){t2[i] <- 2.4*(i-0.5)/n}
f2 <- -t2*sin(-t2^2)
#Scenario 2, k=15, CL=10%, n=300
n <- 300
f_Ridge3 <- rowMeans(f_Ridgex[[18]])
f_Lasso3 <- rowMeans(f_Lassox[[18]])
f_aLasso3 <- rowMeans(f_aLassox[[18]])
f_SCAD3 <- rowMeans(f_SCADx[[18]])
f_MCP3 <- rowMeans(f_MCPx[[18]])
f_Enet3 <- rowMeans(f_Enetx[[18]])
t3 <- 0
for (i in 1:n){t3[i] <- 2.4*(i-0.5)/n}
f3 <- -t3*sin(-t3^2)
#Scenario 2, k=15, CL=30%, n=300
n <- 300
f_Ridge4 <- rowMeans(f_Ridgex[[22]])
f_Lasso4 <- rowMeans(f_Lassox[[22]])
f_aLasso4 <- rowMeans(f_aLassox[[22]])
f_SCAD4 <- rowMeans(f_SCADx[[22]])
f_MCP4 <- rowMeans(f_MCPx[[22]])
f_Enet4 <- rowMeans(f_Enetx[[22]])
t4 <- 0
for (i in 1:n){t4[i] <- 2.4*(i-0.5)/n}
f4 <- -t4*sin(-t4^2)
#------------------------------------------------------------------------------
dat <- simdata(n,0.10,15,1,1)
synz1 <- dat$syn.z
fps1 <- 0
for (i in 1:n){
  if (dat$delta[i]==0){
    fps1[i] <- dat$y[i]-dat$X[i,]%*%dat$beta 
  }
  else{
    fps1[i] <- synz[i]-dat$X[i,]%*%dat$beta
  }
}
fp1 <- f1+rnorm(n)*0.5
df_point1 <- data.frame(t1,fp1)
df_synp1 <- data.frame(t1,fps1)
df_real1 <- data.frame(t1,f1)
df_lasso1 <- data.frame(t1,(f_Lasso1))
df_ridge1 <- data.frame(t1,(f_Ridge1))
df_alasso1 <- data.frame(t1,(f_aLasso1))
df_scad1 <- data.frame(t1,(f_SCAD1))
df_mcp1 <- data.frame(t1,(f_MCP1))
df_enet1 <- data.frame(t1,(f_Enet1))

p1 <- ggplot()+geom_point(data=df_synp1,aes(x=t1,y=fps1/10),col="red",shape=17,alpha=0.4)+geom_point(data=df_point1,aes(x=t1,y=fp1),alpha=0.4)+geom_line(data=df_lasso1,aes(x=t1,y=f_Lasso1-0.2,colour="f(LL)"),lwd=1)+geom_line(data=df_ridge1,aes(x=t1,y=f_Ridge1+0.25,colour="f(RL)"),lwd=1)+geom_line(data=df_alasso1,aes(x=t1,y=f_aLasso1-0.45,colour="f(aLL)"),lwd=1)+geom_line(data=df_scad1,aes(x=t1,y=f_SCAD1,colour="f(SL)"),lwd=1)+geom_line(data=df_mcp1,aes(x=t1,y=f_MCP1,colour="f(MCL)"),lwd=1)+geom_line(data=df_enet1,aes(x=t1,y=f_Enet1,colour="f(ENL)"),lwd=1)+geom_line(data=df_real1,aes(x=t1,y=f1,colour="f(Real)"),lwd=1)+xlab("t")+ylab("f(t) & synthetic data")+ggtitle("(a) Scenario 1, n=50, k=15, CL=10%")+theme(legend.title = element_blank())
p1
#------------------------------------------------------------
dat <- simdata(n,0.30,40,1,1)
synz2 <- dat$syn.z
fps2 <- 0
for (i in 1:n){
  if (dat$delta[i]==0){
    fps2[i] <- dat$y[i]-dat$X[i,]%*%dat$beta 
  }
  else{
    fps2[i] <- synz2[i]-dat$X[i,]%*%dat$beta
  }
}
fp2 <- f2+rnorm(n)*0.5
df_point2 <- data.frame(t2,fp2)
df_synp2 <- data.frame(t2,fps2)
df_real2 <- data.frame(t2,f2)
df_lasso2 <- data.frame(t2,(f_Lasso2))
df_ridge2 <- data.frame(t2,(f_Ridge2))
df_alasso2 <- data.frame(t2,(f_aLasso2))
df_scad2 <- data.frame(t2,(f_SCAD2))
df_mcp2 <- data.frame(t2,(f_MCP2))
df_enet2 <- data.frame(t2,(f_Enet2))

p2 <- ggplot()+geom_point(data=df_synp2,aes(x=t2,y=fps2/4),col="red",shape=17,alpha=0.4)+geom_point(data=df_point2,aes(x=t2,y=fp2),alpha=0.4)+geom_line(data=df_lasso2,aes(x=t2,y=f_Lasso2,colour="f(LL)"),lwd=1)+geom_line(data=df_ridge2,aes(x=t2,y=f_Ridge2,colour="f(RL)"),lwd=1)+geom_line(data=df_alasso2,aes(x=t2,y=f_aLasso2,colour="f(aLL)"),lwd=1)+geom_line(data=df_scad2,aes(x=t2,y=f_SCAD2,colour="f(SL)"),lwd=1)+geom_line(data=df_mcp2,aes(x=t2,y=f_MCP2,colour="f(MCL)"),lwd=1)+geom_line(data=df_enet2,aes(x=t2,y=f_Enet2,colour="f(ENL)"),lwd=1)+geom_line(data=df_real2,aes(x=t2,y=f2,colour="f(Real)"),lwd=1)+xlab("t")+ylab("f(t) & synthetic data")+ggtitle("(b) Scenario 1, n=150, k=40, CL=30%")+theme(legend.title = element_blank())
p2
#------------------------------------------------------------------
#Scenario 2, k=15, CL=10%, n=50
n <- 300
dat <- simdata(n,0.10,15,2,1)
synz3 <- dat$syn.z
fps3 <- 0
for (i in 1:n){
  if (dat$delta[i]==0){
    fps3[i] <- dat$y[i]-dat$X[i,]%*%dat$beta 
  }
  else{
    fps3[i] <- synz3[i]-dat$X[i,]%*%dat$beta
  }
}
fp3 <- f3+rnorm(n)*0.5
df_point3 <- data.frame(t3,fp3)
df_synp3 <- data.frame(t3,fps3)
df_real3 <- data.frame(t3,f3)
df_lasso3 <- data.frame(t3,(f_Lasso3))
df_ridge3 <- data.frame(t3,(f_Ridge3))
df_alasso3 <- data.frame(t3,(f_aLasso3))
df_scad3 <- data.frame(t3,(f_SCAD3))
df_mcp3 <- data.frame(t3,(f_MCP3))
df_enet3 <- data.frame(t3,(f_Enet3))

p3 <- ggplot()+geom_point(data=df_synp3,aes(x=t3,y=fps3),col="red",shape=17,alpha=0.4)+geom_point(data=df_point3,aes(x=t3,y=fp3),alpha=0.4)+geom_line(data=df_lasso3,aes(x=t3,y=f_Lasso3+0.15,colour="f(LL)"),lwd=1)+geom_line(data=df_ridge3,aes(x=t3,y=f_Ridge3,colour="f(RL)"),lwd=1)+geom_line(data=df_alasso3,aes(x=t3,y=f_aLasso3-0.4,colour="f(aLL)"),lwd=1)+geom_line(data=df_scad3,aes(x=t3,y=f_SCAD3-0.15,colour="f(SL)"),lwd=1)+geom_line(data=df_mcp3,aes(x=t3,y=f_MCP3,colour="f(MCL)"),lwd=1)+geom_line(data=df_enet3,aes(x=t3,y=f_Enet3-0.2,colour="f(ENL)"),lwd=1)+geom_line(data=df_real3,aes(x=t3,y=f3,colour="f(Real)"),lwd=1)+xlab("t")+ylab("f(t) & synthetic data")+ggtitle("(c) Scenario 2, n=300, k=15, CL=10%")+theme(legend.title = element_blank())
p3
#-------------------------------------------------------------------------------
#Scenario 2, k=15, CL=30%, n=300
n <- 300
dat <- simdata(n,0.30,15,2,1)
synz4 <- dat$syn.z
fps4 <- 0
for (i in 1:n){
  if (dat$delta[i]==0){
    fps4[i] <- dat$y[i]-dat$X[i,]%*%dat$beta 
  }
  else{
    fps4[i] <- synz4[i]-dat$X[i,]%*%dat$beta
  }
}
fp4 <- f4+rnorm(n)*0.5
df_point4 <- data.frame(t4,fp4)
df_synp4 <- data.frame(t4,fps4)
df_real4 <- data.frame(t4,f4)
df_lasso4 <- data.frame(t4,(f_Lasso4))
df_ridge4 <- data.frame(t4,(f_Ridge4))
df_alasso4 <- data.frame(t4,(f_aLasso4))
df_scad4 <- data.frame(t4,(f_SCAD4))
df_mcp4 <- data.frame(t4,(f_MCP4))
df_enet4 <- data.frame(t4,(f_Enet4))

p4 <- ggplot()+geom_point(data=df_synp4,aes(x=t4,y=fps4/1.2),col="red",shape=17,alpha=0.4)+geom_point(data=df_point4,aes(x=t4,y=fp4),alpha=0.4)+geom_line(data=df_lasso4,aes(x=t4,y=f_Lasso4-0.1,colour="f(LL)"),lwd=1)+geom_line(data=df_ridge4,aes(x=t4,y=f_Ridge4+0.25,colour="f(RL)"),lwd=1)+geom_line(data=df_alasso4,aes(x=t4,y=f_aLasso4-0.35,colour="f(aLL)"),lwd=1)+geom_line(data=df_scad4,aes(x=t4,y=f_SCAD4,colour="f(SL)"),lwd=1)+geom_line(data=df_mcp4,aes(x=t4,y=f_MCP4,colour="f(MCL)"),lwd=1)+geom_line(data=df_enet4,aes(x=t4,y=f_Enet4-0.1,colour="f(ENL)"),lwd=1)+geom_line(data=df_real4,aes(x=t4,y=f4,colour="f(Real)"),lwd=1)+xlab("t")+ylab("f(t) & synthetic data")+ggtitle("(d) Scenario 2, n=400, k=15, CL=30%")+theme(legend.title = element_blank())
p4

grid.arrange(p1,p2,p3,p4,nrow=2)

MSE_all <- matrix(0,24,10)
for (j in 1:24){
  MSE_all[j,] <- colMeans(MSEx[[j]])
}
colnames(MSE_all) <- c("Ridge","Lasso","aLasso","SCAD","MCP","Enet","Scen","n","CL","k")
MSE_all
