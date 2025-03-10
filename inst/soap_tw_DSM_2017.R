#Script soap_tw_DSM_2017.R...Megan C. Ferguson...27 January 2025

  #Notes
  #
  # 1. This script builds the soap film smoother density surface model for 
  #    Eastern Bering Sea belugas in 2017 from Ferguson et al. (2025). This
  #    DSM has the following characteristics:
  #     a. Bivariate (x,y) soap film smoother is isotropic, but respects barriers
  #     b. Tweedie pdf for counts
  #     c. log link
  # 
  # 2. Required input files:
  #    a. Data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata
  #    b. R/mcf_mod_eval_plots.R
  #    c. cpp/soap_tw_DSM.cpp
  #    d. cpp/null_tw_DSM.cpp
  #
  # 3. This script requires the R package TMB. For information on 
  #    installing TMB, see https://github.com/kaskr/adcomp/wiki/Download
  #
  # 4. Figures are output to a folder called "Figures" in the working directory.
  #
  # 5. This script is a clean version of MCF's NSDL17dsm_tmb.R.

    library(mgcv)
    library(TMB)
    library(Matrix)
    library(dsm)
    library(DHARMa)
    library(tweedie)
    library(mgcViz)

    #This script has an option to set k in the argument for the s(bs="so") in the gam.
    #The variable set.k may be set to any value. 
    #The SPDE Matern for this year had 199 random effects (i.e., the
    #number of columns in the A matrix was 199). Setting k=150 (resulting in 199
    #spline coefficients plus one Intercept) is the closest that I can get to putting 
    #the soap model on equal ground for comparison with the SPDE Matern model. 
      set.k <- 150
      
    #Input necessary objects
    
      load("Data/FergusonEtal_20250125_EBS_Beluga_DSM_data.Rdata")
      source("R/mcf_mod_eval_plots.R")
      #CK
        summary(predgrid.strat)
        summary(gam.dat17)
        summary(strat.soap.knots.x)
        class(strat.soap.knots.x) #data.frame
        summary(strat.bnd.list)
        class(strat.bnd.list) #list

    #Create mgcv model
      
      b <- gam(formula = seg.ind ~ s(x, y, bs="so", k=set.k, xt=list(bnd=list(strat.bnd.list))) +
                                             offset(log(a.p)),
                                      family=tw(link="log"), 
                                      method="REML",
                                      data=gam.dat17, 
                                      knots=strat.soap.knots.x)  
      length(b$coefficients) #should be 200
            
    #Evaluate mgcv model fit using 1) function mod.eval.plots from mcf_mod_eval_plots.R,
    #which saves plots to the location specified by the fnam argument to 
    #mod.eval.plots; and 2) DHARMa plots
            
      Fnam <- "Figures/soap_tw_DSM_2017"
            
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

    #Define data object that is given to TMB
                
      #Create the hat matrix in mgcv
              
        Lp <- predict(b, newdata=predgrid.strat, type="lpmatrix")
              
      #Create list of data  
              
        data_tmb = list(y=gam.dat17$seg.ind,      # Response
                        S1=b$smooth[[1]]$S[[1]],  # Penalty matrix for boundary
                        S2=b$smooth[[1]]$S[[2]],  # Penalty matrix for internal
                        X=model.matrix(b)[, -1],  # X*beta is the smoother
                        X_pred=Lp[, -1],          # design matrix for density surface predictions
                        offset=gam.dat17$a.p,     # offset term on the scale of the response
                        offset_pred=predgrid.strat$a.p,  # offset term for density surface predictions
                        n_space_dim=b$smooth[[1]]$null.space.dim)  # Zero for the soap film (but not for other smoother types)

    #Define parameter object given to TMB
                        
      par = list(mu=0,
                 beta=rep(0, length(coef(b))-1),
                 log_lambda=log(b$sp),#c(0, 0),
                 log_phi=0,
                 finv_power=0)

      n.RE <- length(par$beta) #Number of random effects in model
      n.RE
      #CK
        dim(data_tmb$X) #should be 199
        length(par$beta) #should be 199
           
    #Compile and load DLL        
            
      compile("cpp/soap_tw_DSM.cpp")
      dyn.load(dynlib("cpp/soap_tw_DSM"))
            
    #Construct objective functions with derivatives based on compiled C++
    #template       
      M <- MakeADFun(data=data_tmb, parameters=par, random=c("beta"), DLL="soap_tw_DSM")

    #Optimize the model

      Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
      Upper <- 50
      opt.bc <- nlminb( start=M$par, objective=M$fn, gradient=M$gr, lower=Lower, 
                        upper=Upper, control=list(trace=1, eval.max=1000, 
                                                  iter.max=1000))
            
    #Use epsilon detransformation bias correction algorithm
            
      rep <- sdreport( M, par.fixed=opt.bc$par, bias.correct=TRUE)
      #CK
        M$report()
        rep
        summary(rep, "report")

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
        #CK
          summary(pred_tmb) #Predicted plug-in Nhat per cell from TMB
          summary(pred_mgcv) #Predicted plug-in Nhat per cell from mgcv
              
          sum(pred_tmb)  #Predicted plug-in Nhat in study area from TMB
          sum(pred_mgcv) #Predicted plug-in Nhat in study area from mgcv

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
      dyn.unload(dynlib("cpp/soap_tw_DSM"))

    #Steps for computing PDE for the candidate model: 
    # 1. Fit null model; 
    # 2. Compute sum of squared deviances from null model; 
    # 3. Compute sum of residual squared deviances from candidate model;
    # 4. use sum.null.dev.sq and sum.resid.dev.sq to compute PDE.

      #1. Fit null model

        #Compile .cpp for null model then load DLL 
          
          compile("cpp/null_tw_DSM.cpp")
          dyn.load(dynlib("cpp/null_tw_DSM"))

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
              #CK
                null.M$report()
                sum(null.M$report()$devresid^2 )
  
      #Optimize null model
      
        Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
        Upper <- 50
        null.opt.bc <- nlminb( start=null.M$par, objective=null.M$fn, gradient=null.M$gr, 
                               lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, 
                               iter.max=1000))
      #2. Compute null deviance
        sum.null.dev.sq <- sum(null.M$report()$devresid^2 )
        #CK
          sum.null.dev.sq
          
      #3. Recall that the sum of residual squared deviances from candidate model
      #was computed above
        R1 

      #4. Use sum.null.dev.sq and R1 to compute PDE.
        pde <- 1 - R1/sum.null.dev.sq
        #Compare to mgcv's model
          summary(b)
          pde


              
              
              
              
              
