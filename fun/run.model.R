run.model <- function (x, prm, data, LLK=TRUE){
  
  i <- prm$i
  s <- prm$s

  params <- allocate.pars(x,prm)
  
  # Check if parameters are within bounds
  
  tmp <- t(cbind(prm$bds[1:length(x),], (x)))
  tmp <- diff(tmp[c(1,3,2),])

  
  if ((min(tmp)<0)==TRUE){
    llk<- NA
    return(llk)
  }
  
  ########## Set up line of events
  
  # Final conditions, with eventual DST levels
  p4 <- params

  # Health System Restoration
  p3         <-params
  p3$pMDR_rec <-params$pMDR_rec*params$pxpert
  p3$pDx     <-params$pDx*params$pxpert

  
  # During soviet union collapse - escalation in MDR acqu
  p2         <-params
  p2$pMDR_rec <-c(0 , 0)
  p2$pDx     <- params$pDx*0

  #With MDR, but before recent improvements in DST
  p1 <- params
  p1$pMDR_rec <-c(0 , 0.1)
  p1$pDx     <- params$pDx
  
  # Equilibrium: In absence Tb treatments 
  p0 <- params
  p0$beta_mdr <- 0 
  p0$pMDR_acqu <- 0
  p0$pDx      <- 0

  
  
  ######Solve the models 
  
  #Generate a model
  sys.ode.model   <- gen.ode.model(f=system.ode)
  
  
  # Equilibrium Model
  t0 <- 0
  t1 <- 200
  timesODE <- seq(t0,t1,by=1)
  state0 <- rep(0,i$nx)
  seed <- 1e-6
  state0[i$U]        <- (1-seed) 
  state0[i$I$ds$new] <- seed
  sys.df <- sys.ode.model(state=state0,times=timesODE,params=p0)
  next.init<-as.numeric(sys.df[nrow(sys.df),2:ncol(sys.df)])
  soln0<-sys.df[,2:i$nx+1]
  
  
  # With MDR, but before recent improvements in DST
  t0 <- 1970
  t1 <- 1985
  timesODE <- seq(t0,t1,by=1)
  state0<-next.init
  sys.df <- sys.ode.model(state=state0,times=timesODE,params=p1)
  next.init<-as.numeric(sys.df[nrow(sys.df),2:ncol(sys.df)])
  soln1<-sys.df[,2:i$nx+1]
  
  
  # Collapse of USSR
  t0 <- 1985
  t1 <- 1996
  timesODE <- seq(t0,t1,by=1)
  state0<-next.init
  t.interv   <- c(timesODE[2], timesODE[2]+4)
  fx_hndl<-function(t, state, params) scale.ode.model(t, state, params, p1,t.interv)
  sys.scale.model <- gen.ode.model(f=fx_hndl)
  sys.df<- sys.scale.model(state=state0,times=timesODE,params=p2) 
  next.init<-as.numeric(sys.df[nrow(sys.df),2:ncol(sys.df)])
  soln2<-sys.df[,2:i$nx+1]

  # Health System Restoration
  t0 <- 1996
  t1 <- 2012
  timesODE <- seq(t0,t1,by=1)
  state0<-next.init
  t.interv   <- c(timesODE[2], timesODE[2]+5)
  fx_hndl<-function(t, state, params) scale.ode.model(t, state, params, p2,t.interv)
  sys.scale.model <- gen.ode.model(f=fx_hndl)
  sys.df<- sys.scale.model(state=state0,times=timesODE,params=p3) 
  next.init<-as.numeric(sys.df[nrow(sys.df),2:ncol(sys.df)])
  soln3<-sys.df[,2:i$nx+1]
  
  # DST scale up
  t0 <- 2012
  t1 <- 2018
  timesODE <- seq(t0,t1,by=1)
  state0   <-next.init
  t.interv <- c(timesODE[2], timesODE[2]+3)
  fx_hndl  <-function(t, state, params) scale.ode.model(t, state, params, p3,t.interv)
  sys.scale.model <- gen.ode.model(f=fx_hndl)
  sys.df  <- sys.scale.model(state=state0,times=timesODE,params=p4) 
  soln4    <-sys.df[,2:i$nx+1]

  # % --- Get the objectives --------------------------------------------------
  soln<-soln4
  t   <-sys.df[,1]
  sfin<-soln[nrow(soln),]
  gs   <-function(cols) sum(sfin[cols])
  
  # Model output
  N           <-  gs(s$nstates) 
  
  #Prevalence
  prev<- (gs(s$prevalent)/N)*1e5
  
  #Incidence
  tmp_n  <- apply(soln[,i$aux$inc],2,diff)
  tmp    <- 1e5*( tmp_n[nrow(tmp_n),]/N)
  inc_all<- sum(tmp[c(1,2,3,4)])
  inc_mdr<- sum(tmp[c(2,4,5)])
  pmdr   <- 1e2*c(tmp[2]/(tmp[1]+tmp[2]), (tmp[4]+tmp[5])/(tmp[3]+tmp[4]))
 
  
  #Mort
  tmp  <- apply(soln[,i$aux$mort],2,diff)
  tbmort  <- 1e5*(tmp[nrow(tmp),2]/N)
  
  #Notif
  tmp  <- apply(soln[,i$aux$notif],2,diff)
  notif<-1e5*( sum(tmp[nrow(tmp),])/N)
 
  
 
  
  #Model output to fit to data
  model<-c(prev, inc_all,inc_mdr,tbmort, pmdr, notif)
  
  
  
  if(length(data)==0){data<-model}#else{data[c(1:3)]<-data[c(1:3)]/1000 }
  
  llksum<-  sum( (data - model)^2 )
  
  
  # Uncomment for Likelihood method 
  # llk<-c(1:length(model))*0
  # for (l in 1:length(model)){
  #   
  #   llk[l]<-as.numeric(lhd[[l]](model[l]))
  #   if(is.nan( llk[l])){ llk[l]<-0}
  #   
  # }
  # 
  # 
  # llksum<-sum(llk)
  
  # if (is.na(llksum)){
  #   llksum<- NA
  # }
  # 
  
  #  Choos output from this function : If LLK=TRUE only LLK otherwise model output
  
  if (LLK){ 
    
    return(as.numeric(llksum))
    
  }else{
    
    return( list(llk=as.numeric(llksum), soln0=soln0, t0=t0, soln=soln , t=t ,
                 TBepi=model, sys.df=sys.df))
  }
  
  
  
}