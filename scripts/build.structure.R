
# In this script necessary model structure is built and parameters initialized

#__________________________________________________________________________
#  Specifcy model structure and indexing map 
#__________________________________________________________________________

# Model states
gps<-list( strain =c("ds","mdr"),
           hx    =c("new","ret"))

states0 <- c('U')                                            # Unstructured
states1 <- c('Clo','Chi','C','X')                            # structured by strain
states2 <- c('L','I','Dx','E','R','Tx','Tx2')             # structured by strain andretreatment status


i<-list()
s<-list()
groups<-list(states0)
ref <- index.map(groups, i, s, 0)
groups<-list(states1,gps$strain)
ref<-  index.map(groups, ref$i, ref$s, ref$i$nstates)
groups<-list(states2,gps$strain, gps$hx)
ref<-  index.map(groups, ref$i, ref$s, ref$i$nstates)


# Include the auxiliary states

auxnames <- c('inc', 'notif' , 'mort' ,'succ' )
auxinds  <- c(  5      , 2    ,  2    , 2)
ref$i$aux<-list()
lim<-ref$i$nstates
for (ii in 1:length(auxnames)){
  inds <- lim + (1:auxinds[ii])
  ref$i$aux[[auxnames[ii]]] <- inds
  lim <- inds[length(inds)]
}
ref$i$nx <- lim


# useful vectors of indexes
s<-ref$s
i<-ref$i
s$nstates        <- seq(1:i$nstates)
s$infectious     <- c(s$I, s$Dx, s$E, intersect(s$Tx, s$mdr))
s$infectious.ds  <- intersect(c(s$I,  s$Dx, s$E),s$ds)
s$infectious.mdr <- intersect(c(s$I, s$Dx, s$E, s$Tx),s$mdr)
s$prevalent      <- c(s$I, s$Dx, s$E, s$Tx, s$Tx2)
s$TBmort         <- c(s$I, s$Dx, s$E,s$Tx,s$Tx2)
s$sus$new        <- intersect(c(s$L, s$R),s$new)
s$sus$ret        <- c(s$Clo, s$Chi, s$C, intersect(c(s$L, s$R),s$ret))
ref$s<-s
ref$i<-i

#__________________________________________________________________________
#  Model Parameters  
#__________________________________________________________________________

#Label calibration arameters 
xi <-list()

xnames <-c('beta',
           'beta.mdr',
           'careseeking',
           'careseeking2',
           'pDx',
           'MDR_acqu',
           'pxpert', 
           'pTx_init',
           'pTx2_init',
           'Trans')
           
xnums  <- rep(1,length(xnames))
lim <- 0
for (ii in 1:length(xnames)){
  inds <- lim + (1:xnums[ii])
  xi[[xnames[ii]]] <- inds
  lim <- inds[length(inds)]
}

#Assign bounds for sampling
xi$nx <- lim
bds<-matrix(0,length(xnames),2)
bds[xi$beta,]             <- c(1, 18)
bds[xi$beta.mdr,]         <- c(1, 18)
bds[xi$careseeking,]      <- c(0, 12)
bds[xi$careseeking2,]     <- c(9, 14)
bds[xi$pDx,]              <- c(0, 1)  
bds[xi$MDR_acqu,]         <- c(0, 0.1)
bds[xi$pxpert,]           <- c(0, 1)
bds[xi$pTx_init,]         <- c(0.3, 1)
bds[xi$pTx2_init,]        <- c(0.3, 1)
bds[xi$Trans,]            <- c(0, 0.95)   

params<-list(
  
  ####### Calibrated
  beta= 20,
  beta.mdr=4,
  careseeking=2,
  careseeking2=4,
  pDx=0,
  MDR_acqu=0.02,
  pxpert=1, 
  pTx_init=0.9,
  pTx2_init=0.9,
  Trans=0.82,
  
  ####### Fixed
  reactivation=0.001,
  imm=0.5,
  fast=0.14,
  relapse= c(0.032, 0.14, 0.0015, 0.7, 0.02),
  rel_stab=0.5,
  self_cure=1/6,
  rDx=52,
  rTx=2,
  rTx2=0.5,
  pTx_init=0.9,
  pTx2_init=0.9,
  pMDR_rec=c(0.74,0.47),
  cure=c(1,0.3),
  cure2=0.66,
  rdefault=c(0.55,0.66),
  rdefault2=0.44,
  growth=0,
  mu=1/72,
  mu.tb=1/6,

  ####### indexing tools
  pars.n=29,
  lim.pars=-99,
  x.labs=xi,
  bds=bds,
  i=ref$i,
  s=ref$s,
  gps=gps)
