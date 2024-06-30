# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: CoastalWetlands.R
# Author: Lisa Eash
# Date Started: 06/24/2024
# Last Edited: 06/25/2024
# Purpose: AFOLU GHG Calculations for Coastal Wetlands Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: The user must first select the vegetation type (mangrove, tidal marsh, seagrass meadow) and
    #   the subintervention for both the business-as-usual and intervention scenarios. 
    #   Subintervention types include the following:
    #       - Mangrove wood harvest (for harvest or for mangrove management) 
    #       - Rewetting and reestablishment of vegetation
    #       - Draining and vegetation clearing
    #       - Soil excavation and vegetation clearing
    #       - Aquaculture 
    #       - No intervention (only available for business-as-usual scenario)
    #   Other user inputs depend on subintervention category and include the following:
    #       -bio_rem_vol: Removal of mangrove biomass (m^3/ha/year). Only required for "mangrove wood harvest" subintervention
    #       -replanting: Indicates whether reestablishment of vegetation is due to replanting/reseeding or natural recolonization
    #       -salinity: Whether project area is saline (>18 ppt) or freshwater (<18 ppt). Only required for "rewetting and reestablishment of vegetation" subintervention
    #       -fish_production: Annual fish production (kg fish/year). Only required for "aquaculture" subintervention


#### Parameters and Paths ####
  # Calculations requires the following parameters and their associated uncertainty:
      #Biomass C- only account for change in woody biomass (mangroves)
          # Cf: Carbon fraction of mangrove biomass (t C/t dry matter)
          # d: Density of mangrove wood (t/m^3)
          # dom_litter: dead litter C pool for mangrove (t C/ha)
          # dom_wood: dead wood C pool for mangrove (t C/ha)
          # ag_rate: aboveground biomass accumulation rate for mangroves (t dry matter/ha/yr)
          # r: ratio of aboveground to belowground biomass for mangroves
          # ag_biomass: average aboveground biomass for mangroves (t dry matter/ha)
      #SOC
          # EFre: Soil C gains from revegetation (t C/ha/yr; must be from reseeding/planting, not natural recolonization)
          # EFdr: Soil C losses from draining (t C/ha/yr)
          # SObefore: Soil C content in mangroves to calculate losses from soil excavation (t C/ha)
      #N2O
          # EFf: N2O emissions factor for fish production (kg N2O-N per kg fish produced)
      #CH4
          # EFrewet: CH4 emissions factor for soil rewetting (kg CH4/ha/yr)

#### Outputs ####
    # Outputs will be saved to a json file that includes annual CO2, N2O, and CH4 impacts of
    #    the intervention on both a per hectare and total area basis for each subAOI over a 20 year period.
    #    The json file will have the following format:
    # [
    #   {
    #     "AOI": "Sub_AOI-1",
    #     "CO2_thayr": [2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46],
    #     "sdCO2ha": [5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355],
    #     "CO2_tyr": [12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002],
    #     "sdCO2": [28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777],
    #     "N2O_thayr": [2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238],
    #     "sdN2Oha": [3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041],
    #     "N2O_tyr": [11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191],
    #     "sdN2O": [17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206],
    #     "CH4_thayr": [1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491],
    #     "sdCH4ha": [0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849],
    #     "CH4_tyr": [5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457],
    #     "sdCH4": [2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247]
    #   },
    #   {
    #     "AOI": "Sub_AOI-2",
    #     "CO2_thayr": [2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46, 2.46],
    #     "sdCO2ha": [5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355, 5.7355],
    #     "CO2_tyr": [12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002, 12.3002],
    #     "sdCO2": [28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777, 28.6777],
    #     "N2O_thayr": [2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238, 2.3238],
    #     "sdN2Oha": [3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041, 3.4041],
    #     "N2O_tyr": [11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191],
    #     "sdN2O": [17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206],
    #     "CH4_thayr": [1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491],
    #     "sdCH4ha": [0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849],
    #     "CH4_tyr": [5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457],
    #     "sdCH4": [2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247]
    #   }
    # ] 

#### Read in Inputs ####

    # Libraries and wd
    renv::restore()
    library("parallel")
    library("tidyverse")
    library("future.apply")
    library("jsonlite")
    library("rstudioapi")

    current_path = rstudioapi::getActiveDocumentContext()$path 
    setwd(dirname(current_path ))

    # Parameters
    nx=500 #Number of Monte Carlo iterations
    scen=c("business-as-usual","intervention")

    # Extract data from json input
        #Intervention parameters
            j <- jsonlite::fromJSON('coastalwetlandsex.json', flatten=TRUE)
            veg_type<- j$veg_type #vegetation type
            interv_sub_bau<- j$intervention_sub.bau #intervention type for bau scenario
            interv_sub_int<- j$intervention_sub.int #intervention type for intervention scenario
            common<-as.data.frame(j$scenarios$common)
            bau<-as.data.frame(j$scenarios$business_as_usual)
            bau<-bau[,!names(bau)=="aoi_subregions"]
            bau_aoi<-as.data.frame(j$scenarios$business_as_usual$aoi_subregions)
            bau_aoi$aoi_id<-c(1:nrow(bau_aoi))
            bau<-rbind(bau,bau)
            bau<-cbind(bau,bau_aoi)
            bau$scenario<-"business-as-usual"
            int<-as.data.frame(j$scenarios$intervention)
            int<-int[,!names(int)=="aoi_subregions"]
            int_aoi<-as.data.frame(j$scenarios$intervention$aoi_subregions)
            int_aoi$aoi_id<-c(1:nrow(int_aoi))
            int<-rbind(int,int)
            int<-cbind(int,int_aoi)
            int$scenario<-"intervention"
            df<-dplyr::bind_rows(bau, int)
            common<-rbind(common,common,common,common)
            df<-cbind(common,df)
            df <- type.convert(df, as.is = T)
            AOIs<-unique(df$aoi_id)
        #convert uncertainty values to sd
            #Uncertainty presented as 95% CI
            df$cf_sd=(df$cf_uncertainty_high-df$cf_uncertainty_low)/3.92
            df$d_sd=(df$d_uncertainty_high-df$d_uncertainty_low)/3.92
            df$dom_wood_sd=(df$dom_wood_uncertainty_high-df$dom_wood_uncertainty_low)/3.92
            df$dom_litter_sd=(df$dom_litter_uncertainty_high-df$dom_litter_uncertainty_low)/3.92
            df$EFre_sd=(df$EFre_uncertainty_high-df$EFre_uncertainty_low)/3.92
            df$EFdr_sd=(df$EFdr_uncertainty_high-df$EFdr_uncertainty_low)/3.92
            df$EFf_sd=(df$EFf_uncertainty_high-df$EFf_uncertainty_low)/3.92            
            df$EFrewet_sd=(df$EFrewet_uncertainty_high-df$EFrewet_uncertainty_low)/3.92
            df$ag_rate_sd=(df$ag_rate_uncertainty_high-df$ag_rate_uncertainty_low)/3.92
            df$r_sd=(df$r_uncertainty_high-df$r_uncertainty_low)/3.92
            df$ag_biomass_sd=(df$ag_biomass_uncertainty_high-df$ag_biomass_uncertainty_low)/3.92
            df$SObefore_sd=(df$SObefore_uncertainty_high-df$SObefore_uncertainty_low)/3.92
            #Remove other uncertainty columns 
            df<-df[,-grep("uncertainty",names(df))]


####GHG Calculations####
    GHGcalc<-function(i,df,nx #number of MC iterations
    ){
      temp_bau<-df[df$aoi_id==i & df$scenario=="business-as-usual",]
      temp_int<-df[df$aoi_id==i & df$scenario=="intervention",]
      results.unc <- matrix(0, nx, 15)
      results.unc[, 1] <- i
      results.unc[, 2] <-temp_bau$area
      results.unc[, 3] <- c(1:nx)
      colnames(results.unc) <- c('AOI','Area','Rep',
                                 'CO2_y1_bau',
                                 'CO2_y2_bau',
                                 'N2O_bau','CH4_bau',
                                 'CO2_y1_int',
                                 'CO2_y2_int',                            
                                 'N2O_int','CH4_int',
                                 'CO2_y1',
                                 'CO2_y2',
                                 'N2O','CH4')
      for(m in c(1:nx)){
        ####Select parameters####
        cf<-rnorm(1, temp_bau$cf, temp_bau$cf_sd)
        d<-rnorm(1, temp_bau$d, temp_bau$d_sd)
        dom_wood<-rnorm(1, temp_bau$dom_wood, temp_bau$dom_wood_sd)
        dom_litter<-rnorm(1, temp_bau$dom_litter, temp_bau$dom_litter_sd)
        EFre<-rnorm(1, temp_bau$EFre, temp_bau$EFre_sd)
        EFdr<-rnorm(1, temp_bau$EFdr, temp_bau$EFdr_sd)
        EFf<-rnorm(1, temp_bau$EFf, temp_bau$EFf_sd)
        EFrewet<-rnorm(1, temp_bau$EFrewet, temp_bau$EFrewet_sd)
        ag_rate<-rnorm(1, temp_bau$ag_rate, temp_bau$ag_rate_sd)
        r<-rnorm(1, temp_bau$r, temp_bau$r_sd)
        ag_biomass<-rnorm(1, temp_bau$ag_biomass, temp_bau$ag_biomass_sd)
        SObefore<-rnorm(1, temp_bau$SObefore, temp_bau$SObefore_sd)
        for(j in c(1:length(scen))){
          s<-scen[j]
          temp<-df[df$aoi_id==i & df$scenario==s,]
          interv_sub<-ifelse(s=="business-as-usual",interv_sub_bau, interv_sub_int)
          ####Biomass C#### 
              #Woody C gains due to existing or reestablished mangrove
              if(grepl("mangrove wood harvest|reestablishment|no intervention", interv_sub) &&
                 grepl("mangrove",veg_type)){
                Cg=ag_rate*(1+r)*cf
              } else{ Cg=0 }
              #Woody C losses due to harvest or mangrove clearing
              if(grepl("mangrove wood harvest", interv_sub) &&
                 grepl("mangrove",veg_type)){
                Cl=temp$bio_rem_vol*d*cf
              } else{ Cl=0 }
              #Woody C losses due to mangrove clearing - only in first year
              if(grepl("clearing", interv_sub) && #only in first year
                 grepl("mangrove",veg_type)){
                Cly1=ag_biomass*(1+r)*cf
                } else { Cly1=0 }
              #Calculate gain-loss difference for first year and subsequent years
              BioCy1=Cl+Cly1-Cg
              BioCy2=Cl-Cg
              #Dead OM C - reestablishment or clearing of mangroves- just in first year
              if(grepl("clearing", interv_sub) &&
                 grepl("mangrove",veg_type)){
                DOMCy1=dom_litter + dom_wood
              } else { 
                if(grepl("reestablishment", interv_sub) &&
                   grepl("mangrove",veg_type)){
                  DOMCy1=-(dom_litter + dom_wood)
                } else { DOMCy1=0 } }
          ####SOC#### 
              if(grepl("reestablishment", interv_sub) &&
                 temp$replanting=="replanted or reseeded"){
                SOC=EFre
              } else { 
                if(grepl("draining", interv_sub)){
                  SOC=EFdr
                } else { SOC = 0}}
              if(grepl("excavation", interv_sub)){
                  SOCy1=-SObefore
              } else { SOCy1 = 0}
          ####CH4 from rewetting soils#### 
              if(grepl("reestablishment", interv_sub)){
                CH4=EFrewet
              } else { CH4=0 }
          ####N2O from aquaculture####
              if(grepl("aquaculture", interv_sub)){
                N2O=(temp$fish_prod)*EFf*44/28
              } else { N2O=0 }
          #Record in matrix
              results.unc[m,4+(j-1)*4]<-(BioCy1 + DOMCy1 + SOC + SOCy1)*44/12
              results.unc[m,5+(j-1)*4]<-(BioCy2 + SOC)*44/12
              results.unc[m,6+(j-1)*4]<-(N2O)/1000 #kg to t
              results.unc[m,7+(j-1)*4]<-(CH4)/1000 #kg to t
        }
      }
        results.unc[,12]<-results.unc[,8]-results.unc[,4]
        results.unc[,13]<-results.unc[,9]-results.unc[,5]
        results.unc[,14]<-results.unc[,10]-results.unc[,6]
        results.unc[,15]<-results.unc[,11]-results.unc[,7]
        results.unc<-results.unc[,c(1:3,12:15)]
        results.unc<-as.data.frame(results.unc)
        return(results.unc)
    }
      
    #Set up parallel computing
        num_cores <- detectCores()
        cl <- makeCluster(num_cores)
        plan(multisession, workers = num_cores)
        clusterExport(cl, varlist = c("df"))
    #Loop through sub AOIs
        set.seed(1)
        mc<-future_lapply(
          1:length(AOIs), 
          function (x) GHGcalc(x, df, nx),
          future.seed=TRUE
        )
        
    #Extract mean and sd from each subAOI and store in results df 
        resdf <- matrix(0, length(AOIs), 18)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area','CO2_thayr_y1','CO2_thayr_y2',
                             'N2O_thayr','CH4_thayr',
                         'sdCO2ha_y1','sdCO2ha_y2','sdN2Oha','sdCH4ha',
                         'CO2_tyr_y1','CO2_tyr_y2','N2O_tyr','CH4_tyr',
                         'sdCO2_y1','sdCO2_y2','sdN2O','sdCH4')
        for(a in c(1:length(AOIs))){
          resdf[a, 2] <-mean(mc[[a]]$Area)
          resdf[a, 3] <- mean(mc[[a]]$CO2_y1)
          resdf[a, 4] <- mean(mc[[a]]$CO2_y2)
          resdf[a, 5] <- mean(mc[[a]]$N2O)
          resdf[a, 6] <- mean(mc[[a]]$CH4)
          resdf[a, 7] <- sd(mc[[a]]$CO2_y1)
          resdf[a, 8] <- sd(mc[[a]]$CO2_y2)
          resdf[a, 9] <- sd(mc[[a]]$N2O)
          resdf[a, 10] <- sd(mc[[a]]$CH4)
          resdf[a, 11] <- resdf[a,2]*resdf[a,3]
          resdf[a, 12] <- resdf[a,2]*resdf[a,4]
          resdf[a, 13] <- resdf[a,2]*resdf[a,5]
          resdf[a, 14] <- resdf[a,2]*resdf[a,6]
          resdf[a, 15] <- resdf[a,2]*resdf[a,7]
          resdf[a, 16] <- resdf[a,2]*resdf[a,8]
          resdf[a, 17] <- resdf[a,2]*resdf[a,9]
          resdf[a, 18] <- resdf[a,2]*resdf[a,10]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        restib<-data.frame()
        for(r in c(1:nrow(resdf))){
          aoi_tib <-tibble(
            AOI = paste0("Sub_AOI-",r),
            CO2_thayr = list(c(resdf[r,2],rep(resdf[r,3],19))),
            sdCO2ha = list(c(resdf[r,6],rep(resdf[r,7],19))),
            CO2_tyr = list(c(resdf[r,10],rep(resdf[r,11],19))),
            sdCO2 = list(c(resdf[r,14],rep(resdf[r,15],19))),
            N2O_thayr = list(rep(resdf[r,4],20)),
            sdN2Oha = list(rep(resdf[r,8],20)),
            N2O_tyr = list(rep(resdf[r,12],20)),
            sdN2O = list(rep(resdf[r,16],20)),
            CH4_thayr = list(rep(resdf[r,5],20)),
            sdCH4ha  = list(rep(resdf[r,9],20)),
            CH4_tyr = list(rep(resdf[r,14],20)),
            sdCH4  = list(rep(resdf[r,17],20))
          )
          restib<-rbind(restib,aoi_tib)
        }
        json_output<-jsonlite::toJSON(restib,
                                      auto_unbox = T,
                                      pretty = T)
        write(json_output, paste0("output/",interv_sub,".json"))
