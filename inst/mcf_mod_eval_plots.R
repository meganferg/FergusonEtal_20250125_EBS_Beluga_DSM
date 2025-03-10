#Script mcf_mod_eval_plots.R...Megan C. Ferguson...8 May 2024

#NOTES
#1. This script works with quasipoisson, poisson, Tweedie, tw, nb, or negbin models 
#   built using gam() function from package mgcv.
#2. Plots produced:
#   i. Linear predictor vs. randomized quantile residuals
#   ii. QQ plot of deviance residuals
#   iii. Theoretical variance vs. mean, and average squared response residuals 
#        vs. fitted values, within intervals of the predicted values. For 
#        background info about this plot, see:
#          Ver Hoef, Jay M., and Peter L. Boveng. 2007. “QUASI-POISSON VS. 
#          NEGATIVE BINOMIAL REGRESSION: HOW SHOULD WE MODEL OVERDISPERSED 
#          COUNT DATA?” Ecology 88 (11): 2766–72. 
#          https://doi.org/10.1890/07-0043.1.
#3. The arguments are as follows:
#   m: gam model object
#   mod.typ: "nb" or "negbin" = negative binomial, "qp" = quasipoisson, 
#            "tw" or "Tweedie" = Tweedie, "p"=poisson
#   tw.p: Power of the Tweedie distribution, needed only for mod.typ "Tweedie"
#   fnam: output filename

      mod.eval.plots <- function(m, mod.typ, tw.p=NA, fnam){
        
          #1. linear predictor vs. randomized quantile resids
            png(paste(fnam,"_LP_vs_RQR.png",sep=""),bg="white")
              rqgam_check(m)
            dev.off()  
            
          #2. QQ plot of deviance resids  
            png(paste(fnam,"_QQ_dev.png",sep=""),bg="white")
              qq.gam(m, rep = 0, level = 0.9, type = "deviance", rl.col = 2, 
                   rep.col = "gray80")
            dev.off()  
            
          #3. squared response residuals vs. mean
            
            #Compute predictions based on original data
              p <- predict.gam(m, type="response", na.action=na.omit)
            
            #Divide predictions into bins
              phist <- hist(p, breaks=20, plot=FALSE)
              pint <- findInterval(p, phist$breaks)
            
            #Sequence of means for plotting
              mu <- seq(from=range(p)[1], to=range(p)[2], length.out=100)
            
            #Compute squared residuals and average within bins
              #Response residuals
                modl.rr <- residuals(m, type="response")
                modl.rr2 <- (modl.rr)^2
                mean.modl.rr2 <- sapply(1:max(pint), function(i){
                  rr2.i <- modl.rr2[pint == i]
                  mean.rr2.i <- mean(rr2.i)
                  return(mean.rr2.i)
                })
                
            #Compute estimated variances based on mu or p and estimated scale parameter
              if(mod.typ == "nb"){
                v <- mu + (mu^2)/m$family$getTheta(TRUE)
                v.all <- p + (p^2)/m$family$getTheta(TRUE)
              }else if(mod.typ == "negbin"){
                v <- mu + (mu^2)/m$family$getTheta()
                v.all <- p + (p^2)/m$family$getTheta()        
              }else if(mod.typ == "qp" | mod.typ == "p") {
                v <- mu*summary(m)$scale
                v.all <- p*summary(m)$scale
              }else if(mod.typ == "Tweedie") { #V = phi*mu^p
                v <- (summary(m)$scale)*mu^tw.p
                v.all <- (summary(m)$scale)*p^tw.p
              }else{ #"tw"
                v <- (summary(m)$scale)*mu^m$family$getTheta(TRUE)
                v.all <- (summary(m)$scale)*p^m$family$getTheta(TRUE)
              }  
                
            #Plot binned results
              png(paste(fnam,"MuVarAvg_rr.png",sep=""),bg="white")
                plot(phist$mids, mean.modl.rr2, pch=19, xlab="mu",
                     ylab="Mean Squared Response Residuals")
                points(mu, v, pch="*", col="blue")
              dev.off()

      }
