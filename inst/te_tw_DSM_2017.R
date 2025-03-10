#Script te_tw_DSM_2017.R...Megan C. Ferguson...10 March 2025

  #Notes
  #
  # 1. This script builds the tensor product smoother density surface model for 
  #    Eastern Bering Sea belugas in 2017 from Ferguson et al. (2025). This
  #    DSM has the following characteristics:
  #     a. tensor product smooth of x and y, comprising thin-plate smoothing 
  #        splines with shrinkage
  #     b. Tweedie pdf for counts
  #     c. log link
  # 
  # 2. Required input files:
  #    a. data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata
  #    b. inst/mcf_mod_eval_plots.R
  #    c. src/te_tw_DSM.cpp
  #    d. src/null_tw_DSM.cpp
  #
  # 3. This script requires the R package TMB. For information on 
  #    installing TMB, see https://github.com/kaskinst/adcomp/wiki/Download
  #
  # 4. Figures are output to a folder called "Figures" in the working directory.
  #
  # 5. This script is a clean version of MCF's NSDL17dsm_tmb_tw_te.R.

    library(mgcv)
    library(TMB)
    library(Matrix)
    library(dsm)
    library(DHARMa)
    library(tweedie)
    library(mgcViz)

    #   This script has an option to set k in the argument for te( ) in the gam.
    #   The variable set.k may be set to any value. 
    #
    #   The "chosen" SPDE Matern for this year had 199 random effects (i.e., the
    #   number of columns in the A matrix was 199). The total number of coefficients 
    #   in a bivariate te( ) smooth equals k^2, and sqrt(199) = 14.1, so setting k=14
    #   is the closest that I can get to putting the te( ) model on equal ground
    #   for comparison with the SPDE Matern model.      
      set.k <- 14
      
    #Input necessary objects
    
      load("data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata")
      source("inst/mcf_mod_eval_plots.R")

    #Create mgcv model
      
      b <- gam(formula = seg.ind ~ te(x, y, bs="ts", k=set.k) +
                                             offset(log(a.p)),
                                      family=tw(link="log"), 
                                      method="REML",
                                      data=gam.dat17)   
      length(b$coefficients) #should be 196
            
    #Evaluate mgcv model fit using 1) function mod.eval.plots from mcf_mod_eval_plots.R,
    #which saves plots to the location specified by the fnam argument to 
    #mod.eval.plots; and 2) DHARMa plots
            
      Fnam <- "figures/te_tw_DSM_2017"
            
      mod.eval.plots(m=b, mod.typ="tw", tw.p=NA, fnam=Fnam)
            
      #DHARMa function simulateResiduals( )  
            
        b.sim <- simulateResiduals(fittedModel = b, method="PIT", 
                                       n=1000, seed=123)
              
      #Plot: These residuals should be approximately uniform on (0,1) for a well-
      #fitting model.

        png(paste(Fnam, "_DHARMa_mgcv.png",sep=""),bg="white",
                      height=640, width=960, units="px", pointsize=16)
                    plot(b.sim)
                dev.off()  
                
    #Set up only spline structure, following example in pSplines_mcf.R
      gam_setup <- gam(formula = seg.ind ~ te(x, y, bs="ts", k=set.k) +
                                       offset(log(a.p)),
                                family=tw(link="log"), 
                                method="REML",
                                data=gam.dat17,
                                fit = FALSE)

    #Extract penalization matrices.
    #  Because the first (and only) smooth is a bivariate smooth of x and y, 
    #  it has 2 penalty matrices associated with it.
      S_xy_1 <- gam_setup$smooth[[1]]$S[[1]]
      S_xy_2 <- gam_setup$smooth[[1]]$S[[2]]
        
      #Convert to sparse matrices
        S1 <- as(S_xy_1, "CsparseMatrix")
        S2 <- as(S_xy_2, "CsparseMatrix")

      #Extract number of parameters needed to define the smooth
        n.beta <- b$smooth[[1]]$last.para - b$smooth[[1]]$first.para + 1 ## number of params
        n.beta #195
            
      #Extract the number of spline penalization coeffs
        n.lambda <- length(b$smooth[[1]]$S)
        n.lambda #2

    #Define data object that is given to TMB
                
      #Create the hat matrix in mgcv
              
        Lp <- predict(b, newdata=predgrid.strat, type="lpmatrix")
              
      #Create list of data  
              
        data_tmb <- list(y=gam.dat17$seg.ind,  # Response
                        X = gam_setup$X[,-1],  # Design matrix, without intercept
                        X_pred=Lp[, -1],       # design matrix for density surface predictions
                        offset = gam.dat17$a.p,     # offset term on the scale of the response
                        offset_pred=predgrid.strat$a.p,  # offset term for density surface predictions
                        S1 = S1,      # sparse penalty mtx for first term
                        S2 = S2)      # sparse penalty mtx for second term

      #Define parameter object given to TMB
                        
        par = list(
                  beta0 = 0,  # Intercept
                  beta = rep(0,n.beta),  #Spline coefficients. There is a single model matrix for 
                                         #a tensor product smooth. 
                  log_lambda = rep(0,n.lambda), #Log spline penalization coefficients.
                                                #One penalty for each marginal basis in
                                                #a tensor product smooth.
                  log_phi=0,
                  finv_power=0
                  )

        n.RE <- length(par$beta) #Number of random effects in model
        n.RE

    #Compile and load DLL        
            
      compile("src/te_tw_DSM.cpp")
      dyn.load(dynlib("src/te_tw_DSM"))
            
    #Construct objective functions with derivatives based on compiled C++
    #template       
      M <- MakeADFun(data = data_tmb, 
                     parameters = par, 
                     random="beta",
                     DLL = "te_tw_DSM")

    #Optimize the model

      Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
      Upper <- 50
      opt.bc <- nlminb( start=M$par, objective=M$fn, gradient=M$gr, lower=Lower, 
                        upper=Upper, control=list(trace=1, eval.max=1000, 
                                                  iter.max=1000))
            
    #Use epsilon detransformation bias correction algorithm
            
      rep <- sdreport( M, par.fixed=opt.bc$par, bias.correct=TRUE)

    #Compare cell-wise predictions (plug-in estimator) from mgcv and TMB. 

      #TMB
              
        idx <- which(names(summary(rep)[,1]) == "log_lambda" |
                     names(summary(rep)[,1]) == "log_phi" |
                     names(summary(rep)[,1]) == "finv_power" |   
                     names(summary(rep)[,1]) == "Nhat" |
                     names(summary(rep)[,1]) == "Nhat_pred")
        pred_tmb <- predgrid.strat$a.p*exp(Lp%*%summary(rep)[-idx, 1]) 

      #mgcv point estimate of Nhat per cell via Gaussian approximation 
      #(plug-in estimate)  
                
        pred_mgcv <- predgrid.strat$a.p*exp(Lp%*%coef(b)) #log link, with offset

    #Evaluate TMB model fit           
              
      #Extract estimates of phi and power for Tweedie
          
        #TMB        
          idx <- which(names(summary(rep)[,1]) == "log_phi")   
          tmb.phi <- exp(summary(rep)[idx,1])  
            
          idx <- which(names(summary(rep)[,1]) == "finv_power")  
          finv.power <- summary(rep)[idx,1]
          tmb.power <- 1.0 + (exp(finv.power) / (1.0 + exp(finv.power)))
            
        #mgcv  
          summary(b)$scale #phi from mgcv
          tmb.phi #should be similar to mgcv
            
          b$family$getTheta(TRUE) #power from mgcv
          tmb.power #should be nearly identical to mgcv
     
      #Compute fitted values
          
        #Find betas
          
          idx <- which(names(summary(rep)[,1]) == "log_lambda" |
                             names(summary(rep)[,1]) == "log_phi" |
                             names(summary(rep)[,1]) == "finv_power" |   
                             names(summary(rep)[,1]) == "Nhat" |
                             names(summary(rep)[,1]) == "Nhat_pred")
          
        #Create the hat (design) matrix in mgcv and generate tmb fits
          Lp.fit.tmb <- cbind(rep(1, nrow(data_tmb$X)), data_tmb$X)
          fit.tmb <- data_tmb$offset*exp(Lp.fit.tmb%*%summary(rep)[-idx, 1]) 
    
      #Compute and plot residuals
        
        #Simulate response (replicated data sets) from the predictive distribution 
        #of data conditional on estimated fixed and random effects. 
          n.sims <- 1000 #250 is the default for DHARMa function simulateResiduals( )
              
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
            
        #Plot tmb model residuals
          png(paste(Fnam, "_DHARMa_TMB.png",sep=""),bg="white",
                height=640, width=960, units="px", pointsize=16)
              plot(DHARMa.res)
          dev.off() 
              
    #Compute sum of residual squared deviances from candidate model. This
    #will be used later in the script to calculate the percent deviance explained.
      R1 = sum(M$report()$devresid^2 )
              
    #Unload DLL for candidate DSM     
      dyn.unload(dynlib("src/te_tw_DSM"))

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
    
          null.data_tmb = list(y=gam.dat17$seg.ind,  # Response
                              offset=gam.dat17$a.p)  # offset term on the scale of the response
  
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


              
              
              
              
              
