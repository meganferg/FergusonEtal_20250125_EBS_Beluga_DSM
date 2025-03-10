#Script EBS_beluga_Nht_2022.R...Megan C. Ferguson...10 March 2025

  #Notes
  #
  # 1. This script derives an estimate of EBS beluga abundance from 2022
  #    aerial line-transect survey data using the Horvitz-Thompson estimator
  #    described in Ferguson et al. (2025).
  #
  # 2. Required input files:
  #    a. data/FergusonEtal_20250125_EBS_Beluga_Nht_data.Rdata
  #
  # 6. This script is based on MCF's NSDL2022abund.R.

    library(Distance) 

    #Input stuff
      load("data/FergusonEtal_20250125_EBS_Beluga_Nht_data.Rdata")

    #Compute Nhat for full 2022 study area  

        Dl.strat.D.4.6 <- dht(model=Dl.trnc5pct.4.hr.iBeauf.strat6, 
                             region.table=Region.Table.6,
                             sample.table=smpl.strat.table.4.6,
                             obs.table=Dl.x95.strat.4.6,
                             se=TRUE)

    #Compute Nhat for portion of 2022 study area included in 2017 estimate 
      Dl.L4.D.4 <- dht(model=Dl.trnc5pct.4.hr.iBeauf.strat, 
                               region.table=Region.Table.L4,
                               sample.table=smpl.strat.table.L4,
                               obs.table=L4.obs,
                               se=TRUE)

  #Use delta method to incorporate cv.p0.1 into total CV
    cv.Nhat.4.6 <- sqrt( (cv.p0.1^2) + (Dl.strat.D.4.6$individuals$N[7,4])^2 )    
    cv.Nhat.L4 <- sqrt( (cv.p0.1^2) + (Dl.L4.D.4$individuals$N[5,4])^2 )    
