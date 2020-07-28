# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# R Code for TB epidemic in Euro Region and covid impact
# Author: Juan F Vesga (2020)
#
# Simulation functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model ODE as required by deSolve
system.ode <- function(t,state,params) {
  with(as.list(c(state, params)), {
    
    invec <- as.numeric(state[1:i$nstates])
    N     <- sum(invec)
    morts <- sum(invec*mu)  + sum(invec[s$TBmort]*mu.tb)  
    births<- morts*(growth==0) + growth
    foi<- c(beta*sum(invec[s$infectious.ds])/N , beta.mdr*sum(invec[s$infectious.mdr])/N)      
    
    # if (t>2){
    #   browser()
    #   }
    # 
    dx<-state*0
    
    dx[i$U] <- births - invec[i$U]*sum(foi) - invec[i$U]*mu
    
    
    for (is in 1:length(gps$strain)){
      strain <- gps$strain[is]
      ismdr   <- strcmp(strain, 'mdr')
      
      Clo  <- i$Clo[[strain]]
      Chi  <- i$Chi[[strain]]
      C    <- i$C[[strain]]
      X    <- i$X[[strain]]
      
      
      for (ih in 1:length(gps$hx)){
        hx <- gps$hx[ih]
        isret   <- strcmp(hx, 'ret')
        
        L   <- i$L[[strain]][[hx]]
        I   <- i$I[[strain]][[hx]]
        Dx  <- i$Dx[[strain]][[hx]]
        E   <- i$E[[strain]][[hx]]
        R   <- i$R[[strain]][[hx]]
        Tx  <- i$Tx[[strain]][[hx]]
        Tx2 <- i$Tx2[[strain]][[hx]]  
        
        
        dx[L] <- invec[i$U]*foi[is]*(1-fast)*(1-isret) + 
          sum(invec[s$sus[[hx]]])*foi[is]*(1-fast)*imm -
          invec[L]*reactivation - 
          invec[L]*(sum(foi))*imm - 
          invec[L]*mu
        
        
        dx[I] <- invec[L]*reactivation +
          invec[i$U]*foi[is]*fast*(1-isret) +
          sum(invec[s$sus[[hx]]])*foi[is]*fast*imm +
          isret*(invec[Clo]*relapse[1] + invec[Chi]*relapse[2] +
                   invec[C]*relapse[3]+invec[X]*relapse[4]) + 
          invec[R]*relapse[5]  - invec[I]*(self_cure+careseeking+mu+mu.tb)
        
        
        # --- Diagnosis 
        pFLinit <- pDx*pTx_init*(1 - ismdr*pMDR_rec[ih])
        pSLinit <- pDx*pTx2_init*ismdr*pMDR_rec[ih]
        p_ltfu  <- 1- (pFLinit+pSLinit)
        
        dx[Dx] <- invec[I]*careseeking + invec[E]*careseeking2 -
          invec[Dx]*rDx -
          invec[Dx]*(self_cure+mu+mu.tb)
        

        # --- FL Treatment
        pFLcure <- cure[is]
        pSLtran <- Trans*ismdr
        
        
        dx[Tx] <- invec[Dx]*rDx*pFLinit + 
          invec[i$Tx$ds[[hx]]]*MDR_acqu*ismdr - 
          invec[Tx]*(rTx+rdefault[is]+MDR_acqu*(1-ismdr)) -
          invec[Tx]*(mu+mu.tb)
        
        dx[Tx2] <- invec[Dx]*rDx*pSLinit + 
          invec[i$Tx$mdr$new]*rTx*(1-pFLcure)*pSLtran*isret +
          invec[i$Tx$mdr$ret]*rTx*(1-pFLcure)*pSLtran*isret - 
          invec[Tx2]*(rTx2+rdefault2) -
          invec[Tx2]*(mu+mu.tb)
        
        dx[E] <- invec[Dx]*rDx*p_ltfu + 
          invec[i$Tx$mdr$new]*rTx*(1-pFLcure)*(1-pSLtran)*isret +
          invec[i$Tx$mdr$ret]*rTx*(1-pFLcure)*(1-pSLtran)*isret + 
          invec[Tx]*rdefault[is] +
          invec[Tx2]*rTx2*(1-cure2)*isret + 
          invec[Tx2]*rdefault2*isret -
          invec[E]*careseeking2 - 
          invec[E]*(self_cure+mu+mu.tb)
        
        dx[R] <- sum(invec[c(I,E,Dx)])*self_cure - 
          invec[R]*(sum(foi))*imm - 
          invec[R]*(relapse[5]+mu)
        
        #outcomes
        dx[i$aux$notif[1]]<-dx[i$aux$notif[1]] + invec[Dx]*rDx*pFLinit
        dx[i$aux$notif[2]]<-dx[i$aux$notif[2]] + invec[Dx]*rDx*pSLinit
        
        
      }
      
      # --- Relapse from therapeutic cure compartments
      dx[Clo]<- sum(invec[as.numeric(i$Tx$ds)]*rTx*cure[1])*(1-ismdr) +
        sum(invec[as.numeric(i$Tx2$mdr)]*rTx2*cure2)*(ismdr) -
        invec[Clo]*(sum(foi))*imm -
        invec[Clo]*(relapse[1]+rel_stab+mu)
      
      dx[Chi]<- -invec[Chi]*(relapse[2]+rel_stab+mu) - 
        invec[Chi]*(sum(foi))*imm 
      
      dx[C]  <- sum(invec[c(Clo,Chi)])*rel_stab - 
        invec[C]*(sum(foi))*imm -
        invec[C]*(relapse[3]+mu)
      
      dx[X]  <- sum(invec[as.numeric(i$Tx$mdr)]*rTx*cure[2])*(ismdr) - 
        invec[X]*relapse[4]
      
      dx[i$aux$succ[1]] <-  dx[i$aux$succ[1]] + sum(invec[as.numeric(i$Tx$ds)]*rTx*cure[1])*(1-ismdr)
      dx[i$aux$succ[2]] <- dx[i$aux$succ[2]] + sum(invec[as.numeric(i$Tx2$mdr)]*rTx2*cure2)*(ismdr)
      
    }
    
    
    #outcomes
    dx[i$aux$inc[1]] <- invec[i$U]*foi[1]*fast + 
      sum(invec[s$sus[["new"]]])*foi[1]*fast*imm + 
      sum(invec[i$R$ds$new])*relapse[5] +
      sum(invec[i$L$ds$new])*reactivation
    
    dx[i$aux$inc[2]] <- invec[i$U]*foi[2]*fast + 
      sum(invec[s$sus[["new"]]])*foi[2]*fast*imm + 
      invec[i$R$mdr$new]*relapse[5] +
      invec[i$L$mdr$new]*reactivation
    
    dx[i$aux$inc[3]] <- sum(invec[s$sus[["new"]]])*foi[1]*fast*imm +   
      invec[i$Clo$ds]*relapse[1] + 
      invec[i$Ch$ds]*relapse[2] +
      invec[i$C$ds]*relapse[3]+
      invec[i$X$ds]*relapse[4] +
      invec[i$R$ds$ret]*relapse[5] +
      invec[i$L$ds$ret]*reactivation
    
    dx[i$aux$inc[4]] <- sum(invec[s$sus[["new"]]])*foi[2]*fast*imm +   
      invec[i$Clo$mdr]*relapse[1] + 
      invec[i$Ch$mdr]*relapse[2] +
      invec[i$C$mdr]*relapse[3]+
      invec[i$X$mdr]*relapse[4] +
      invec[i$R$mdr$ret]*relapse[5] +
      invec[i$L$mdr$ret]*reactivation
    
    dx[i$aux$inc[5]] <- sum(invec[c(i$Tx$ds$new, i$Tx$ds$ret)])*MDR_acqu
    
    dx[i$aux$mort[1]]<- sum(invec*mu)   
    dx[i$aux$mort[2]]<- sum(invec[s$TBmort]*mu.tb) 
    
    
    
    
    
    return(list(dx))
  })
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Scale parameters model
scale.ode.model <- function(t,state, p1, p0,t.interv) {
  
  scale <- min((t-t.interv[1])/(t.interv[2]-t.interv[1]),1) 
  
  if (scale<0) 
  {
    scale<-0
  }
  # lim<-which(p0 %in% -99)-1
  pars_old<-p0[1:p0$pars.n]
  pars_new<-p1[1:p0$pars.n]
  
  pars_now<- Map(function(num0,num1,scale) num0+scale*(num1-num0),pars_old,pars_new,scale)
  

  pars_scaled<-p1
  pars_scaled[1:p0$pars.n]<-pars_now
  
  return(system.ode(t,state,pars_scaled))
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate a model (wrapper)
gen.ode.model <- function(f) {
  fret <- function(state,times,params) {
    ode.out <- lsoda(y = state, times = times, func = f, parms = params)
    return(as.data.frame(ode.out))
  }
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# auxiliary functions
gen.plot.ode <- function(traj.df) {
  traj.m = melt(traj.df, id.vars='time')
  p <- ggplot(traj.m, aes(x = time, y = value, color = variable)) +
    theme_bw() + geom_line()
  return(p)
}