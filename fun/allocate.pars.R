allocate.pars<-function(x,prm){
  
  xi<- prm$x.labs
  
  prm$beta          <-  x[xi$beta]
  prm$beta.mdr      <-  x[xi$beta.mdr]
  prm$careseeking   <-  x[xi$careseeking]
  prm$careseeking2  <-  x[xi$careseeking2]
  prm$pDx           <-  x[xi$pDx]
  prm$MDR_acqu      <-  x[xi$MDR_acqu]
  prm$pxpert        <-  x[xi$pxpert]
  prm$pTx_init      <-  x[xi$pTx_init]
  prm$pTx2_init     <-  x[xi$pTx2_init]
  prm$Trans         <-  x[xi$Trans]
  
  return(prm)
}


