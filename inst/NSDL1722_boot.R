#Script NSDL1722_boot.R...Megan C. Ferguson...24 August 2024

# 1. This script is based on Conn's run_jobss_boot.R
#
# 2. In the TMB model bootstrap, dyn.load only needs to occur once, outside the 
#    bootstrap. No need to dyn.unload at the end of each bootstrap iteration. 

  library(TMB)
  library(mgcv)
  library(Matrix)
  library(mrds)
  
  setwd("C:\\Users\\megan.ferguson\\OneDrive - Biodiversity Research Institute\\Documents\\Belugas\\NSDLabund\\Analysis\\NSDLabund_BRI")
  
  #Specify path to cpp files
    cpp.path <- "C:\\Users\\megan.ferguson\\cpp_temp\\"
  
  #Set number of bootstrap iterations
    nboot <- 500
    
    #nboot <- 10
  
  #Detection Function
  
    #Input fitted detection function, segment data, prediction grid
      load("Output_DSM//Ferguson_NSDL1722dsm_stuff.Rdata")
      ddf.obj <- Dl.trnc5pct.hr.iBeauf.Turb

    #Extract the MLEs and covariance matrix from fitted detection function model
      ddf.beta <- Dl.trnc5pct.hr.iBeauf.Turb$par
      ddf.Sigma <- solve(ddf.obj$hessian)
      #CK
        summary(ddf.obj)
        
        ddf.beta #These values should equal the parameter estimates in the 
                  #summary call above. First value is for shape coeff.
        
        sqrt(diag(ddf.Sigma)) #These values should equal the SEs for the 
                              #parameter estimates in the summary(ddf) line above.
                              #First value is for shape coeff.
        
    #Implement ddf bootstrap
        
      set.seed(7777)
        
      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      #Create matrices to hold bootstrapped values of average.p. Rows
      #correspond to segments and columns correspond to bootstrap iterations.
        bs.p17 <- matrix(data=rep(0, nboot*nrow(seg.dat17.x)),
                             nrow=nrow(seg.dat17.x))
        
        bs.p22 <- matrix(data=rep(0, nboot*nrow(seg.dat22)),
                             nrow=nrow(seg.dat22))
        #CK
          dim(seg.dat17.x)
          dim(bs.p17)
          
          dim(seg.dat22)
          dim(bs.p22)
        
      for(ddf.boot in 1:nboot){
        
        #ddf.boot <- 1

        cat(paste0("ddf.boot = ",ddf.boot,"\n"))
        
        #Sample ddf betas, assuming a multivariate normal distribution
          ddf.beta.i <- mgcv::rmvn(1, ddf.beta, ddf.Sigma)
          #CK
            if(DEBUGG){
              print(ddf.beta)
              print(ddf.beta.i)
            }
          
        #generate new bootstrap detection probabilities
          bs.ddf <- ddf.obj
          bs.ddf$par <- ddf.beta.i
          
          bs.p17[,ddf.boot] <- predict(bs.ddf, seg.dat17.x)$fitted 
          bs.p22[,ddf.boot] <- predict(bs.ddf, seg.dat22)$fitted 
          
      }
      #CK
        dim(seg.dat17.x)
        dim(bs.p17)
        
        summary(bs.p17[1,])
        predict(ddf.obj, seg.dat17.x[1,])$fitted #Should be similar to above
          
        dim(seg.dat22)
        dim(bs.p22)
        
        summary(bs.p22[nrow(seg.dat22),])
        predict(ddf.obj, seg.dat22[nrow(seg.dat22),])$fitted #Should be similar to above

      

    
  ###########      
  #TMB model#
  ###########      
        
    M.name <- "NSDL17dsm_tmb_tw_soap_xy_sf"
    DLL.name <- "NSDL17dsm_tmb_PDE"
    Rdat.name <- "NSDL17dsm_tmb_tw_soap_xy_w_spde_knots_ns20000_k105"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
        
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data_tmb 
          bs.dat$offset <- gam.dat17$area * bs.p17[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random=c("beta"), DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)

        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[2,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[2,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 
      
  ###########    
  #TMB model#
  ###########    
        
    M.name <- "NSDL17dsm_tmb_mgcv_spde_tw_mesh6_sf"
    DLL.name <- "NSDL17dsm_tmb_spde_tw_PDE"
    Rdat.name <- "NSDL17dsm_tmb_mgcv_spde_tw_mesh6"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
      
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data 
          bs.dat$offset <- seg.dat17.x$area * bs.p17[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random="x", DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[3,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[3,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 
      
  ###########    
  #TMB model#
  ###########    
    
    M.name <- "NSDL17dsm_tmb_tw_te_sf"
    DLL.name <- "NSDL17dsm_tmb_tw_te_PDE"
    Rdat.name <- "NSDL17dsm_tmb_tw_te_ns20000_k14"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
      
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data 
          bs.dat$offset <- gam.dat17$area * bs.p17[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random="beta", DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[2,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[2,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
  ###########    
  #TMB model#
  ###########    
    
    M.name <- "NSDL17dsm_tmb_tw_xy_sf"
    DLL.name <- "NSDL17dsm_tmb_tw_x_y_PDE"
    Rdat.name <- "NSDL17dsm_tmb_tw_xy_ns20000_k200"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
      
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data 
          bs.dat$offset <- gam.dat17$area * bs.p17[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random="beta", DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[2,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[2,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
  ###########    
  #TMB model#
  ###########    
        
    M.name <- "NSDL22dsm_tmb_mgcv_spde_tw_noscmn_mesh5_sf"
    DLL.name <- "NSDL17dsm_tmb_spde_tw_PDE"
    Rdat.name <- "NSDL22dsm_tmb_mgcv_spde_tw_noscmn_mesh5"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
      
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data 
          bs.dat$offset <- seg.dat22$area * bs.p22[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random="x", DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[3,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[3,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 
      
  ###########    
  #TMB model#
  ###########    
    
    M.name <- "NSDL22dsm_tmb_mgcv_spde_tw_bnd_noscmn_mesh2_sf"
    DLL.name <- "NSDL17dsm_tmb_spde_tw_bnd_PDE"
    Rdat.name <- "NSDL22dsm_tmb_spde_tw_bnd_noscmn_mesh2"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
      
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data 
          bs.dat$offset <- seg.dat22$area * bs.p22[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random="x", DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[3,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[3,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 
      
  ###########      
  #TMB model#
  ###########      
        
    M.name <- "NSDL22dsm_tmb_tw_soap_xy_w_spde_knots_noscmn_sf"
    DLL.name <- "NSDL17dsm_tmb_PDE"
    Rdat.name <- "NSDL22dsm_tmb_tw_soap_xy_w_spde_knots_noscmn_ns20000_k165"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
        
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data_tmb 
          bs.dat$offset <- gam.dat22$area * bs.p22[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random=c("beta"), DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[2,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[2,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 

  ###########  
  #TMB model#
  ###########    

    M.name <- "NSDL22dsm_tmb_tw_te_noscmn_sf"
    DLL.name <- "NSDL17dsm_tmb_tw_te_PDE"
    Rdat.name <- "NSDL22dsm_tmb_tw_te_noscmn_ns20000_k18"
        
    #Load compiled cpp code

      dyn.load(dynlib(paste(cpp.path, DLL.name, sep=""))) 
    
    #Input necessary R objects from previously saved workspace.

      load(paste("Output_DSM//",Rdat.name,".Rdata",sep="")) 
      
    #Bootstrap  
      
      #Prep dataframe to save bootstrap Nhat and se values
      
        #nboot <- 2 #For debugging
      
        set.seed(20230929)
        
        Nhat <- rep(0,nboot)
        se <- rep(0,nboot)
        conv <- rep(NA,nboot)
        boot.df <- cbind.data.frame(Nhat, se, conv) #bias-corrected Nhat and SE(Nhat) from each BS iter

      DEBUGG <- FALSE
      #DEBUGG <- TRUE
      
      for(iboot in 1:nboot){
        
        #iboot <- 1

        cat(paste0(M.name, " iboot = ",iboot,"\n"))
        
        #Generate bootstrap data list using bootstrapped detection probabilities
        
            ###############################################################
            #Need to tweak variable names for each M.name because I wasn't# 
            #consistent in naming in the original model scripts           #
            ###############################################################        
        
          bs.dat <- data 
          bs.dat$offset <- gam.dat22$area * bs.p22[,iboot] #This is the offset for the GAM.
                                                #It DOES account for p0.1 & pavail.
          
        #Re-fit and optimize TMB model
          bs.M <- MakeADFun(data=bs.dat, parameters=par, random="beta", DLL=DLL.name)
          Lower <- -50  #trying to prevent -Inf,Inf bounds resulting in nlminb failure (NaN gradient)
          Upper <- 50
          bs.opt.bc <- nlminb( start=bs.M$par, objective=bs.M$fn, gradient=bs.M$gr, lower=Lower, upper=Upper, control=list(trace=1, eval.max=1000, iter.max=1000))
          bs.rep <- sdreport( bs.M, par.fixed=bs.opt.bc$par, bias.correct=TRUE)
         
        #Save bs bias-corrected Nhat and se(Nhat), and convergence message 
          bs.rep.mtx <- summary(bs.rep, "report")
          boot.df$Nhat[iboot] <- bs.rep.mtx[2,3] #Might need to tweak indices depending on M.name
          boot.df$se[iboot] <- bs.rep.mtx[2,4]   #Might need to tweak indices depending on M.name
          boot.df$conv[iboot] <- bs.opt.bc$convergence
          #CK
          if(DEBUGG){
            print(bs.rep.mtx)
            print(boot.df)
          }
          
      }
      
    #Save results 
      save(boot.df, file=paste("Output_DSM//",M.name,"_BS.Rdata",sep=""))
      rm(boot.df)
      dyn.unload(dynlib(paste(cpp.path, DLL.name, sep=""))) 

