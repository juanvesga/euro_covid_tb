rm(list = ls())

#__________________________________________________________________________
#  Load libraries and necesary scripts 
#__________________________________________________________________________
library(deSolve)
library(ggplot2)
library(reshape2)
# library(Matrix)
library(readxl)
library(pracma)
# library(neldermead)
library(profvis)
# library(Rcpp)
# library(abind)
# library(openxlsx)

source("fun/model.functions.R")
source("fun/index.map.R")
source("fun/run.model.R")
source("fun/allocate.pars.R")
source("scripts/build.structure.R")

# source("fun/get_objective.R")
# source("fun/goveqs_basis.R")
# source("fun/scale_up.R")

# source("fun/make_model.R")
# source("fun/return_output.R")
# source("fun/get_input.R")
# source("fun/sim_pack_intervention.R")
# source("fun/Calibration.R")
# source("fun/Simulate.R")
# sourceCpp("fun/compute_dx.cpp")
# sourceCpp("fun/scale_matrix.cpp")



#__________________________________________________________________________
#  run a model
#__________________________________________________________________________



# 
# 
# t0 <- 0
# t1 <- 200
# timesODE <- seq(t0,t1,by=1)
# 
# state0 <- rep(0,ref$i$nx)
# seed <- 1e-6
# state0[ref$i$U] <- (1-seed)
# state0[ref$i$E$ds$new] <- seed
# sys.ode.model <- gen.ode.model(f=system.ode)
# sysODEdf <- sys.ode.model(state=state0,times=timesODE,params=params)
# # 
# soln<-sysODEdf[,s$nstates+1]
# N=rowSums(soln)
# plot(N)
# plotODE <- gen.plot.ode(sysODEdf[,seq(0,37)+1])
# print(plotODE)
# 
# plotODE <- gen.plot.ode(sysODEdf[,c(0,as.numeric(i$I$ds$new))+1])
# print(plotODE)
# 
# plotODE <- gen.plot.ode(sysODEdf[,c(0,as.numeric(i$aux$inc[1]))+1])
# print(plotODE)
# 
# plot( diff(sysODEdf[,as.numeric(i$aux$inc[1])+1])*1e5,type = 'l')
# plot( diff(sysODEdf[,as.numeric(i$aux$inc[2])+1])*1e5,type = 'l')

# 
file <-"data/data.xlsx"
# 
# # Load Data : This is required for calibration and simulation steps
fulldata      <-as.data.frame(read_excel(file))
# 
# ##-----------------------------------------------------------------
# ## 1) Run calibrations
# ##-----------------------------------------------------------------

# Select state
 
state<- "Mol"

# Get inputs

tmp<- fulldata[fulldata$State== state, ]
target.data<-tmp[,-1]

x0<-c(9.89881704 ,4.50534668, 1.26884163, 9.06362617, 0.72692015, 0.07684589, 0.77456720 ,0.74819106, 0.77927747, 0.58638353)

# run.model(x0, params, target.data)
# 
# library(profvis)
# 
# 
# profvis({
#  run<-run.model(x0, params, target.data,LLK=FALSE)
# })
#  
#  
#  tmp_n  <- apply(run$soln[,i$aux$inc],2,diff)
#  plot( rowSums(tmp_n)*1e5,type = 'l')
#  
#  browser()
#  
obj_spx   <-  function(x)
  run.model(x, params,target.data)

avgpar<-rowSums(params$bds)/2
ratio<-avgpar[1]/avgpar
ratio<-ratio*200

x1 <- optim(x0,obj_spx,
            control = list(maxit = 300,
                           trace=TRUE,
                           reltol=1e-2,
                           abstol=0.1,
                           parscale = ratio))
x  <- x1$par




# 
# # Call calibration function (C)
# # Note: If a prevoius calibration (Y0) exists function C will retrieve it and 
# # use it as starting parameters
# C<-Calibration(X_Y0, state, district)
# 
# 
# # Save parameters 
# Y1<-data.frame(t(C$x))
# wb<-loadWorkbook(file)
# writeData(wb, "Y", Y1, startCol = 3, startRow = 2, colNames = FALSE)
# 
# # Save workbook
# saveWorkbook(wb, file, overwrite = TRUE)
# 
# 
# 
# 
# ##-----------------------------------------------------------------
# ##- 2) Run interventions
# ##-----------------------------------------------------------------
# # Load previously saved values (Y) to initialize State variables
# Y<-as.data.frame(read_excel(file, sheet = "Y"))
# tmp<- Y[Y$state== state & Y$district == district,  ]
# 
# Y1 <- as.numeric(tmp[c(3:length(tmp))])
# 
# # Get dto inputs
# tmp<- fulldata[fulldata$state== state & fulldata$district == district,  ]
# X<-tmp[,-c(1,2)]
# 
# # Call intervention function "Simulate"
# sims<-
#   Simulate (
#     X, #  District base characteristics
#     Y1,   # set of best fit parameters
#     
#     # ALL Values belwo correspont to I
#     0.95,      # Target for first-line treatment success
#     
#     0.75,      # Target for second-line treatment success
#     
#     0.5,       # 1-Proportion of diagnoses being done by CB-NAAT vs smear (see change in order vs.) 
#     0.95,      # Proportion of diagnoses having DST result within two weeks of diagnosis
#     
#     0.95,      # Target proportion of private providers to be engaged
#     0.8,       # Proportion of diagnoses using Xpert in private engaged
#     0.85,      # Target treatment outcome in privatre engaged
#     
#     1,          #Population  (< 5 only or all household contacts):  0.1 for childen, 1 for al population
#     4,         # Number of contacts envisaged to be screened per index case
#     0.6,       # Anticipated regimen completion rates
#     
#     0.2,       # Proportion of risk group to be screened per year
#     0,         # Screening algorithm (verbally, or also with X-ray):  Xray =1 ; Verbally =0;
#     0,         # Confirmatory test (smear or Xpert): Xpert =1 ; Smear =0;
#     0.8,       # Proportion of presumptives successfully linked to microbiological testing’
#     0.8        # Proportion of diagnosed cases successfully initiating treatment 
#   )
# 
# # Save results (Z)
# tmp<-c(sims$inc_base,
#        sims$mort_base,
#        sims$inc_itv,
#        sims$mort_itv,
#        sims$ic_fl,
#        sims$ic_sl,
#        sims$ic_sm,
#        sims$ic_xp,
#        sims$ic_xr,
#        sims$icr_all,
#        sims$inc_av_n,
#        sims$mort_av_n)
# 
# Z<-data.frame(t(tmp))
# wb<-loadWorkbook(file)
# 
# id<-which(fulldata$state== state & fulldata$district == district)
# 
# writeData(wb, "Z", Z, startCol = 3, startRow = id+1, colNames = FALSE)
# saveWorkbook(wb, file, overwrite = TRUE)
