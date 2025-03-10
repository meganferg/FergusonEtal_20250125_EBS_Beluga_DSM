#Script xy_tw_DSM_2017.R...Megan C. Ferguson...10 March 2025

  #Notes
  #
  # 1. This script builds the density surface model for Eastern Bering Sea belugas 
  #    in 2017 using the SPDE framework from Ferguson et al. (2025). The DSM
  #    is a purely spatial model. The number of belugas sighted on a segment  
  #    with midpoint located at x is assumed to be Tweedie distributed as: 
  #         y_i ~ Tweedie(eta_i)
  #    where 
  #         ln(eta_i) = X*beta + delta(x_i) + log(offset_i)
  #    Here X is a design matrix for the fixed effects, beta is a vector with the 
  #    covariate regression coefficients, delta is a latent spatial field with Matern 
  #    covariance structure without barriers, and the offset equals 2*L_i*ESW_i.
  #
  # 2. This script requires the mgcv SPDE functions in DLM's mgcv_spde_smooth.R.
  #    On 5.8.24, MCF revised mgcv_spde_smooth.R because it was not correctly 
  #    identifying inla.mesh class objects; I saved my edited version to 
  #    mgcv_spde_smooth_mcf.R. For further information about DLM's script,
  #    see Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic 
  #    Partial Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020). 
  #    https://doi.org/10.1007/s13253-019-00377-z
  #
  # 3. Required input files:
  #    a. data/Ferguson_NSDL17_20240508_spde.Rdata
  #    b. inst/mcf_mod_eval_plots.R
  #    c. inst/mgcv_spde_smooth_mcf.R
  #    d. src/spde_tw_DSM.cpp
  #    e. src/null_tw_DSM.cpp
  #
  # 4. This script requires the R package TMB. For information on 
  #    installing TMB, see https://github.com/kaskinst/adcomp/wiki/Download
  #
  # 5. Figures are output to a folder called "Figures" in the working directory.
  #
  # 6. This script is based on MCF's NSDL17dsm_tmb_spde_tw.R and 
  #    NSDL17dsm_tw_spde.R.

    library(mgcv)
    library(TMB)
    library(Matrix)
    library(dsm)
    library(DHARMa)
    library(tweedie)
    library(mgcViz)
    library(INLA)
    library(fields)

    #Designate path and basic filename for outputting figures
      Fnam <- "figures/spde_tw_DSM_2017"
    
    #Input necessary objects
    
      load("data/Ferguson_NSDL17_20240508_spde.Rdata")
      source("inst/mcf_mod_eval_plots.R")
      source("inst/mgcv_spde_smooth_mcf.R")

    #Define mesh and components representing the  precision matrix

      #Use all segment midpoints    
          
        loc = cbind(seg.dat$x, seg.dat$y)
            
        #From the inla.nonconvex.hull helpfile: 
        #  convex: The desired extension radius. Also determines the smallest 
        #         allowed convex curvature radius. Negative values are interpreted 
        #         as fractions of the approximate initial set diameter.
        boundary = INLA::inla.nonconvex.hull(loc)
        
        #From the inla.mesh.2d helpfile:
        #  max.edge: The largest allowed triangle edge length. One or two values. 
        #  cutoff: The minimum allowed distance between points. Point at most as far apart as 
        #          this are replaced by a single vertex prior to the mesh refinement step.
        #
        #See also https://rpubs.com/jafet089/886687 for a vignette with advice and 
        #guidelines on building the mesh.
        
          boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35) 
          max.edg <- 60 #km

          mesh <- INLA::inla.mesh.2d(
            loc=loc,
            boundary=list(boundary,boundary2),
            max.edge=c(1,2)*max.edg,
            cutoff = max.edg/5
          )

    #Continue with spde model specification      

      A = inla.spde.make.A(mesh,loc) #spatial interpolation mtx for observations
        
      spde = inla.spde2.matern(mesh, alpha=2)
        
      spdeMatrices = spde$param.inla[c("M0","M1","M2")]
        
    #Create spatial interpolation matrix for predicting to the grid
        
      pred_loc <- cbind(predgrid.strat$x, predgrid.strat$y)  
      A_pred <- inla.spde.make.A(mesh, pred_loc) #spatial interpolation mtx for predictions  

    #Define the design matrix for the fixed effects
          
      X <- model.matrix( ~ 1,
                         data = seg.dat) #No covariates

    #Define the design matrix for prediction
          
      X_pred <- model.matrix( ~ 1,
                         data = predgrid.strat) #No covariates

    #Define the data and parameters given to TMB

      data_tmb = list(y = seg.dat$seg.ind, #Response
                  A = A,            #Spatial interpolation matrix for observations
                  A_pred = A_pred,  #Spatial interpolation mtx for predictions
                  spdeMatrices = spdeMatrices, #SPDE matrices
                  X = as.matrix(X), #Design matrix for observations
                  X_pred = as.matrix(X_pred), #Design mtx for prediction
                  offset = seg.dat$a.p, #Offset = 2 * L * ESW * p.avail * p0
                  offset_pred=predgrid.strat$off.set  # offset term for density surface predictions
      )
              
      par = list(beta = rep(0,dim(X)[2]),
                 log_tau =1,
                 log_kappa = -3,
                 x = rep(0,mesh$n),
                 log_phi=0,
                 finv_power=0
                )
              
    #Compile and load c++ code
      
      compile("src/spde_tw_DSM.cpp")
      dyn.load(dynlib("src/spde_tw_DSM"))
            
    #Build the model
      
      M <- MakeADFun(data_tmb, par, random="x", DLL="spde_tw_DSM")

    #Optimize the model

      Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
      Upper <- 50
      opt.bc <- nlminb( start=M$par, objective=M$fn, gradient=M$gr, lower=Lower, 
                        upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
            
    #Use epsilon detransformation bias correction algorithm. 
            
      rep <- sdreport( M, par.fixed=opt.bc$par, bias.correct=TRUE)

    #Make cellwise predictions from SPDE model
              
      rangeIndex <- which(row.names(summary(rep,"report"))=="range")
      range.tmb <- summary(rep,"report")[rangeIndex,]
      range.tmb

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

      Nhat.pred.i <- predgrid.strat$off.set*exp(mu.pred)

    #Evaluate TMB model fit using DHARMa          
          
      #Extract estimates of phi and power for Tweedie 
        idx <- which(names(summary(rep)[,1]) == "log_phi")   
        tmb.phi <- exp(summary(rep)[idx,1])  
          
        idx <- which(names(summary(rep)[,1]) == "finv_power")  
        finv.power <- summary(rep)[idx,1]
        tmb.power <- 1.0 + (exp(finv.power) / (1.0 + exp(finv.power)))

      #Compute fitted values

        delta <- as.matrix((A%*%x.est)/tau.est)

        X.beta <- X*beta.est

        mu.tmb <- X.beta + delta
          
        fit.tmb <- data_tmb$offset*exp(mu.tmb)

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
                                         observedResponse = data_tmb$y,
                                         fittedPredictedResponse = fit.tmb[,1],
                                         integerResponse = FALSE,
                                         method="PIT",
                                         seed=123)
            
        #Plot tmb models
          png(paste(Fnam,"_DHARMa_TMB.png",sep=""),bg="white",
                    height=640, width=960, units="px", pointsize=16)
                  plot(DHARMa.res)
          dev.off()
              
    #Build SPDE model in mgcv using DLM's functions
        
      ## fit model using gam( ) function
        # disable scale penalty as otherwise smoothing parameters
        # are not directly interpretable
          M.gam <- gam(seg.ind ~ s(x, y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)) +
                               offset(log(a.p)),
                     data = seg.dat,
                     control = gam.control(scalePenalty = FALSE),
                     family=tw(link="log"), 
                     method = "REML")

        ## get hyperparameter estimates
          tau.gam <- M.gam$sp[1]
          kappa.gam <- M.gam$sp[2]
          # compute correlation range (rho) and marginal variance (sigma)
          range.gam <- sqrt(8)/kappa.gam #From Skaug script spdeBarrier.cpp
          alpha <- 2
          nu <- alpha - 1/2
          
        # predict mgcv model to grid
          
          #Option 1: Use the method in DLM's script campylobacterosis.R
        
            pred.df <- cbind.data.frame("x"=predgrid.strat$x, "y"=predgrid.strat$y)
            
            # get design matrix
              Xp <- PredictMat(M.gam$smooth[[1]], data = pred.df)
              
            # add in intercept
              Xp <- cbind(1, Xp)
              
            # compute posterior mean
              predmu <- predgrid.strat$a.p*exp(Xp %*% coef(M.gam))     
            
          #Option 2:   
            
            # Create the hat (design) matrix in mgcv
              Lp <- predict(M.gam, newdata=predgrid.strat, type="lpmatrix")
            
            #Apply offset and predict to grid      
                pred_mgcv <- predgrid.strat$a.p*exp(Lp%*%coef(M.gam))
                mgcv.GA.Nhat_pred <- sum(pred_mgcv)
                
            #Compare results from Option 1 and Option 2    
                
              length(predmu)
              length(pred_mgcv)
                
              summary(predmu - pred_mgcv)
              sum(predmu - pred_mgcv) #They're identical
                
              summary(predmu) #And they're not all zero!
              summary(pred_mgcv)
              
              summary(Nhat.pred.i - pred_mgcv) #TMB - mgcv
              sum(Nhat.pred.i) #predicted area-integrated Nhat from TMB
              sum(pred_mgcv)   #predicted area-integrated Nhat from mgcv

      ## Try fitting model using dsm( ) function
        # disable scale penalty as otherwise smoothing parameters
        # are not directly interpretable
        M.dsm <- dsm(count~s(x, y, bs="spde", k = mesh$n, xt = list(mesh = mesh)),
                                      family=tw(link="log"), method="REML",
                                      control = gam.control(scalePenalty = FALSE,
                                                            keepData = TRUE),
                                      ddf.obj = Dl.trnc5pct.hr.iBeauf.Turb, 
                                      segment.data = seg.dat,
                                      observation.data = obs.dat)

        #Can predict dsm to grid two ways
            
          pred1_dsm <- predict(M.dsm, newdata=predgrid.strat, type="response")
          
          Lp.dsm <- predict(M.dsm, newdata=predgrid.strat, type="lpmatrix")
          pred2_dsm <- predgrid.strat$a.p*exp(Lp.dsm%*%coef(M.dsm)) 
          #CK to confirm that pred_mgcv equals area-integrated predictions from dsm
            summary(pred1_dsm - pred_mgcv[,1]) #Near zero
            summary(pred2_dsm[,1] - pred1_dsm) #Near zero
            
            sum(pred1_dsm - pred_mgcv[,1])     #Near zero
            sum(pred2_dsm[,1] - pred1_dsm)     #Near zero

        ## get hyperparameter estimates
          tau.dsm <- M.dsm$sp[1]
          kappa.dsm <- M.dsm$sp[2]
          # compute correlation range (rho) and marginal variance (sigma)
          range.dsm <- sqrt(8)/kappa.dsm #From Skaug script spdeBarrier.cpp
          alpha <- 2
          nu <- alpha - 1/2
          
            cat("tmb:\n")
            cat("kappa=", kappa.est, "\n")
            cat("tau=", tau.est, "\n")
            cat("range=", range.tmb[1], "\n")
          
            cat("mgcv:\n")
            cat("kappa=", kappa.gam, "\n")
            cat("tau=", tau.gam, "\n")
            cat("range=", range.gam, "\n")
          
            cat("dsm:\n")
            cat("kappa=", kappa.dsm, "\n")
            cat("tau=", tau.dsm, "\n")
            cat("range=", range.dsm, "\n")
            
    #Evaluate mgcv model fit using 1) function mod.eval.plots from mcf_mod_eval_plots.R,
    #which saves plots to the location specified by the fnam argument to 
    #mod.eval.plots; and 2) DHARMa plots
            
      mod.eval.plots(m=M.gam, mod.typ="tw", tw.p=NA, fnam=Fnam)
            
      #DHARMa function simulateResiduals( )  
            
        M.gam.sim <- simulateResiduals(fittedModel = M.gam, method="PIT", 
                                        n=1000, seed=123)
              
        #Plot: These residuals should be approximately uniform on (0,1) for a well-
        #fitting model.

          png(paste(Fnam, "_DHARMa_mgcv.png",sep=""),bg="white",
                height=640, width=960, units="px", pointsize=16)
              plot(M.gam.sim)
          dev.off()  
                
    #Compute sum of residual squared deviances from candidate model. This
    #will be used later in the script to calculate the percent deviance explained.
      R1 = sum(M$report()$devresid^2 )

    #Unload DLL for candidate DSM     
      dyn.unload(dynlib("src/spde_tw_DSM"))

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
    
          null.data_tmb = list(y=seg.dat$seg.ind,  # Response
                              offset=seg.dat$a.p)  # offset term on the scale of the response
  
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


              
              
              
              
              
