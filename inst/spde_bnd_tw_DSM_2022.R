#Script spde_bnd_tw_DSM_2022.R...Megan C. Ferguson...10 March 2025

  #1. This script creates a DSM using the spde framework. The DSM is a purely 
  #   spatial model. The number of belugas sighted on a segment with midpoint 
  #   located at x is assumed to be Tweedie distributed as: 
  #         y_i ~ Tweedie(eta_i)
  #   where 
  #         ln(eta_i) = X*beta + delta(x_i) + log(offset_i)
  #   Here X is a design matrix for the fixed effects, beta is a vector with the 
  #   covariate regression coefficients, delta is a latent spatial field with Matern 
  #   covariance structure with barriers, and the offset equals 2*L_i*ESW_i. 
  #
  #2. This script loads spde_bnd_tw_DSM (created using spde_bnd_tw_DSM.cpp) as a 
  #   dynamically linked library. In order for spde_bnd_tw_DSM.cpp to compile, 
  #   it needs access to the file "barrierTMB.hpp", which I copied from 
  #   https://github.com/skaug/tmb-case-studies/tree/master/spdeBarrier 
  #
  # 3. Required input files:
  #    a. data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata
  #    b. src/spde_bnd_tw_DSM.cpp
  #    c. src/null_tw_DSM.cpp
  #    d. src/barrierTMB.hpp
  #
  # 4. This script requires the R package TMB. For information on 
  #    installing TMB, see https://github.com/kaskr/adcomp/wiki/Download
  #
  # 5. Figures are output to a folder called "figures" in the working directory.
  #
  # 6. This script is based on MCF's NSDL22dsm_tmb_spde_tw_bnd_noscmn.R.

    library(mgcv)
    library(TMB)
    library(Matrix)
    library(dsm)
    library(DHARMa)
    library(tweedie)
    library(INLA)
    library(mrds)
    library(sf)
    library(sp)
    
    #Designate path and basic filename for outputting figures
      Fnam <- "figures/spde_bnd_tw_DSM_2022"

    #Define projection
      nsdl.proj.km <- CRS('+proj=eqdc +lat_1=62.5d
                                  +lat_2=64.5d
                                  +lat_0=63.5d
                                  +lon_0=-164.0d
                                  +x_0=0
                                  +y_0=0
                                  +datum=WGS84
                                  +units=km') #units are in km
      
          
    #Input necessary objects
      load("data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata")

      #Convert units for x and y in gam.dat22 and predgrid from m to km. 
      #The distances and areas in gam.dat22 and predgrid are already in km.
        
        gam.dat22$x.m <- gam.dat22$x
        gam.dat22$x <- gam.dat22$x.m/1000
        
        gam.dat22$y.m <- gam.dat22$y
        gam.dat22$y <- gam.dat22$y.m/1000
        
        predgrid$x.m <- predgrid$x
        predgrid$x <- predgrid$x.m/1000
        
        predgrid$y.m <- predgrid$y
        predgrid$y <- predgrid$y.m/1000
        
        predgrid.sf <- st_as_sf(predgrid, coords=c("x","y"))
        st_crs(predgrid.sf) <- nsdl.proj.km

    #Spatial interpolation matrix      
      
      loc <- cbind(gam.dat22$x, gam.dat22$y)
        
      A = inla.spde.make.A(bnd.mesh,loc)
      
    #Construct barrier.triangles.
      
      tl = length(bnd.mesh$graph$tv[,1]) #extract the number of triangles in the mesh
      posTri = matrix(0, tl, 2) #matrix of zeros with tl rows and 2 columns
      for(t in 1:tl){
        temp = bnd.mesh$loc[bnd.mesh$graph$tv[t, ], ] #Extract coords of the 3 vertices of triange t
        posTri[t,] = colMeans(temp)[c(1,2)] #Compute the center of triangle t
      }
      posTri = SpatialPoints(posTri, proj4string=nsdl.proj.km) #Convert center of all triangles to SpatialPoints
      
      bnd22.noscmn.km.SP <- spTransform(bnd22.noscmn.SP, nsdl.proj.km)
      normal = unlist(over(bnd22.noscmn.km.SP, posTri, returnList=T)) #triangles located in WATER in study area
      barrier.triangles = setdiff(1:tl, normal) #triangles located either on land or not in study area

    #Construct spde object and structure needed in the spde-procedure
    #  inla.barrier.fem() is an internal method producing the Finite Element matrices  
        
      spde = inla.spde2.matern(bnd.mesh, alpha=2)
      spdeMatrices = spde$param.inla[c("M0","M1","M2")] #Matrices needed in Standard spde-procedure
      fem = INLA:::inla.barrier.fem(mesh = bnd.mesh, barrier.triangles = barrier.triangles)
      spdeMatricesBarrier = list(C0 = fem$C[[1]],C1 =fem$C[[2]] ,D0 = fem$D[[1]],D1 = fem$D[[2]],I = fem$I )

    #Create spatial interpolation matrix for predicting to the grid
        
      pred_loc <- cbind(predgrid$x, predgrid$y)  
      A_pred <- inla.spde.make.A(bnd.mesh, pred_loc) #spatial interpolation mtx for predictions  

    #Define the design matrix for the fixed effects
          
      X <- model.matrix( ~ 1,
                         data = gam.dat22) #No covariates

    #Define the design matrix for prediction
          
      X_pred <- model.matrix( ~ 1,
                         data = predgrid) #No covariates

    #Define the data and parameters given to TMB

      data = list(y = gam.dat22$seg.ind, #Response
                  A = A,            #Spatial interpolation matrix
                  A_pred = A_pred,  #Spatial interpolation mtx for predictions
                  spdeMatrices = spdeMatrices, #SPDE matrices
                  spdeMatricesBarrier = spdeMatricesBarrier, #needed to use the barrier option in TMB
                  barrier = 1,
                  c = c(1,0.2), #Values taken from Skaug's spdeBarrier.R
                                #Second element in c provides the scaling factor of the range on land
                  X = as.matrix(X), #Design matrix
                  X_pred = as.matrix(X_pred), #Design mtx for prediction
                  offset = gam.dat22$a.p, #Offset = 2 * L * ESW * p.avail * p0
                  offset_pred = predgrid$off.set  # offset term for density surface predictions
      )

      par = list(beta = rep(0,dim(X)[2]), #vector of covariate regression coeffs
                 log_tau =1,
                 log_kappa = -3,
                 x = rep(0,bnd.mesh$n), #value of the random field at the mesh nodes
                 log_phi=0,
                 finv_power=0
                )
      
      n.RE <- length(par$x)

    #Compile and load c++ code

      compile("src/spde_bnd_tw_DSM.cpp")
      dyn.load(dynlib("src/spde_bnd_tw_DSM"))
            
    #Build the model 
      
      M <- MakeADFun(data, par, random="x", DLL="spde_bnd_tw_DSM")

    #Optimize     
                  
      Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
      Upper <- 50
      opt.bc <- nlminb( start=M$par, objective=M$fn, gradient=M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))

    #Use epsilon detransformation bias correction algorithm. 
  
      rep <- sdreport( M, par.fixed=opt.bc$par, bias.correct=TRUE)

    #Make cellwise predictions from SPDE model
        
      fieldIndex <- which(row.names(summary(rep))=="x")
      x.est <- summary(rep)[fieldIndex,1]
      
      log_tau_Idx <- which(row.names(summary(rep))=="log_tau") 
      tau.est <- exp(summary(rep)[log_tau_Idx,1]) 
      
      log_kappa_Idx <- which(row.names(summary(rep))=="log_kappa") 
      kappa.est <- exp(summary(rep)[log_kappa_Idx,1]) 
      
      betaIndex <- which(row.names(summary(rep))=="beta")
      beta.est <- summary(rep)[betaIndex,1]
      
      Nhat_pred_Idx <- which(row.names(summary(rep))=="Nhat_pred")
      Nhat.pred.est <- summary(rep)[Nhat_pred_Idx,1]

      delta.pred <- as.matrix((A_pred%*%x.est)/tau.est)
        class(delta.pred)
        dim(delta.pred)
        
      X.beta <- X_pred*beta.est
        class(X.beta)
        dim(X.beta)
      
      mu.pred <- X.beta + delta.pred

      Nhat.pred.i <- predgrid$off.set*exp(mu.pred)

    #Evaluate TMB model fit using DHARMa          
          
      #Extract estimates of phi and power for Tweedie 
        idx <- which(names(summary(rep)[,1]) == "log_phi")   
        tmb.phi <- exp(summary(rep)[idx,1])  
        tmb.phi
          
        idx <- which(names(summary(rep)[,1]) == "finv_power")  
        finv.power <- summary(rep)[idx,1]
        tmb.power <- 1.0 + (exp(finv.power) / (1.0 + exp(finv.power)))
        tmb.power

      #Compute fitted values

        delta <- as.matrix((A%*%x.est)/tau.est)

        X.beta <- X*beta.est

        mu.tmb <- X.beta + delta
          
        fit.tmb <- data$offset*exp(mu.tmb)

      #Compute and plot residuals using DHARMa
          
        #Simulate response (replicated datasets) from the predictive distribution 
        #of data conditional on estimated fixed and random effects. 
          n.sims <- 500 #250 is the default for DHARMa function simulateResiduals( )
              
          sims <- t(sapply(1:length(fit.tmb), function(j){
                      y.j <- rtweedie(n = n.sims,
                                      mu = fit.tmb[j], 
                                      phi = tmb.phi, 
                                      power = tmb.power)
                      return(y.j)
                     }))

        #Create DHARMa object
          DHARMa.res <- createDHARMa(simulatedResponse = sims, #simulated response values
                                         observedResponse = data$y,
                                         fittedPredictedResponse = fit.tmb[,1],
                                         integerResponse = FALSE,
                                         method="PIT",
                                         seed=123)
            
        #Plot tmb models
          png(paste(Fnam,"_DHARMa_TMB.png",sep=""),bg="white",
                    height=640, width=960, units="px", pointsize=16)
                  plot(DHARMa.res)
          dev.off()
        
    #Compute sum of residual squared deviances from candidate model. This
    #will be used later in the script to calculate the percent deviance explained.
      R1 = sum(M$report()$devresid^2 )

    #Unload DLL for candidate DSM     
      dyn.unload(dynlib("src/spde_bnd_tw_DSM"))

    #Steps for computing PDE for the candidate model: 
    # 1. Fit null model; 
    # 2. Compute sum of squared deviances from null model; 
    # 3. Compute sum of residual squared deviances from candidate model;
    # 4. use sum.null.dev.sq and sum.resid.dev.sq to compute PDE.

      #1. Fit null model

        #Compile .cpp for null model then load DLL 
          
          compile("src/null_tw_DSM.cpp")
          dyn.load(dynlib("src/null_tw_DSM"))

        #Specify data and parameters
    
          null.data_tmb = list(y=gam.dat22$seg.ind,  # Response
                              offset=gam.dat22$a.p)  # offset term on the scale of the response
  
          null.par = list(mu=0,
                         log_phi=0,
                         finv_power=0)
          
        #Construct objective functions with derivatives based on compiled C++
        #template
        
          null.M <- MakeADFun(data=null.data_tmb, 
                         parameters=null.par, 
                         DLL="null_tw_DSM")

      #Optimize null model
      
        Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
        Upper <- 50
        null.opt.bc <- nlminb( start=null.M$par, objective=null.M$fn, gradient=null.M$gr, 
                               lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, 
                               iter.max=1000))
      #2. Compute null deviance
        sum.null.dev.sq <- sum(null.M$report()$devresid^2 )

      #3. Recall that the sum of residual squared deviances from candidate model
      #was computed above
        R1 

      #4. Use sum.null.dev.sq and R1 to compute PDE.
        pde <- 1 - R1/sum.null.dev.sq

        dyn.unload(dynlib("src/null_tw_DSM"))

    #Now derive DSM predictions for the portion of the study area that was
    #surveyed in 2017
      
      #Set up prediction grid: extract only the cells in the 2017 study area,
      #which are identified in sf22.in17
      
        predgrid22.in17 <- predgrid[which(sf22.in17 == 1),]

      #Create spatial interpolation matrix for predicting to the grid
          
        pred_loc <- cbind(predgrid22.in17$x, predgrid22.in17$y)  
        A_pred <- inla.spde.make.A(bnd.mesh, pred_loc) #spatial interpolation mtx for predictions  

      #Define the design matrix for prediction
            
        X_pred <- model.matrix( ~ 1,
                           data = predgrid22.in17) #No covariates

      #Define the data and parameters given to TMB
  
        data = list(y = gam.dat22$seg.ind, #Response
                    A = A,            #Spatial interpolation matrix
                    A_pred = A_pred,  #Spatial interpolation mtx for predictions
                    spdeMatrices = spdeMatrices, #SPDE matrices
                    spdeMatricesBarrier = spdeMatricesBarrier, #needed to use the barrier option in TMB
                    barrier = 1,
                    c = c(1,0.2), #Values taken from Skaug's spdeBarrier.R
                                  #Secound element in c provides the scaling factor of the range on land
                    X = as.matrix(X), #Design matrix
                    X_pred = as.matrix(X_pred), #Design mtx for prediction
                    offset = gam.dat22$a.p, #Offset = 2 * L * ESW * p.avail * p0
                    offset_pred=predgrid22.in17$off.set  # offset term for density surface predictions
        )
  
        par = list(beta = rep(0,dim(X)[2]),
                   log_tau =1,
                   log_kappa = -3,
                   x = rep(0,bnd.mesh$n),
                   log_phi=0,
                   finv_power=0
                  )
        
      #Load compiled c++ code
  
        dyn.load(dynlib("src/spde_bnd_tw_DSM"))
              
      #Estimate the model and extract results
        
        M22.in17 <- MakeADFun(data, par, random="x", DLL="spde_bnd_tw_DSM")
  
        #Bias correction
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          opt.bc <- nlminb( start=M22.in17$par, objective=M22.in17$fn, gradient=M22.in17$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          
        #Apply epsilon detransformation bias correction factor  
          rep <- sdreport( M22.in17, par.fixed=opt.bc$par, bias.correct=TRUE)

                
                    
                    
                    
                    
                    
   