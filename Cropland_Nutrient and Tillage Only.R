# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: Cropland.R
# Author: Lisa Eash
# Date Started: 06/11/2024
# Last Edited: 06/14/2024
# Purpose: AFOLU GHG Calculations for Annual and Perennial Cropland Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: User inputs depend on subintervention category (tillage reduction or
    #   nutrient management). For tillage reduction, only tillage type is required. For
    #   nutrient management, user must provide organic matter management (low, med,
    #   high, high with organic amendments) and synthetic N application rate and N%. If  
    #   the addition of organic amendments is part of the improved nutrient management,
    #   they may also enter organic amendment application rate and N%.

#### Parameters and Paths ####
    # Calculations requires the following parameters and their associated uncertainty:
      #SOC
          # SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
          # FLU: Land use emissions factor for climate zone and soil type
          # FMG: Tillage factor for both business as usual and intervention scenario (if 
          #    tillage is not a subintervention then FMG = 1)
          # FI: C input factor for both business as usual and intervention scenario (if 
          #    nutrient management is not a subintervention then FI = 1)
      #N2O
          # FSN: Amount of synthetic N applied annually (kg N/ha/yr)
          # FON: Amount of organic N applied annually (kg N/ha/yr)
          # EFDir: Emissions factor for direct N2O emissions from N inputs
          # FRAC_GASF: Fraction of synthetic fertilizer that volatilizes
          # FRAC_GASM: Fraction of organic N that volatilizes
          # EF_vol: Emissions factor for indirect N2O emissions from volatilization
          # FRAC_LEACH: Fraction of N inputs that is lost to leaching
          # EF_leach: Emissions factor for indirect N2O emissions from volatilization

#### Outputs ####
    # Outputs will be saved to a json file that includes annual CO2 and N2O impacts of
    #    the intervention on both a per hectare and total area basis for each subAOI

#### Read in Inputs ####

# Libraries
library(jsonlite)
library(tidyverse)
library(parallel)
library(future.apply)
library(rstudioapi)
# it be good to track the specific versions of these library you're using, 
# either manually in the readme or using something like ren if you're familiar 
# I couldn't run this at first because I didn't see the plyr dependency for instance


    # Extract data from json input
        #trying to still set the working directory but more generically, do you think this is okay
        current_path = rstudioapi::getActiveDocumentContext()$path 
        setwd(dirname(current_path ))
        #GEE output
        
            # I don't think you should need this file, because this information 
            # is already in croplandex.json in the aoi_subregions. 
            # The other input file you're getting according to  
            # https://mathematicampr.atlassian.net/wiki/spaces/AFOLU/pages/3212444049/Code+requirements
            # is land cover information which isn't in croplandex.json 
            gee<-jsonlite::fromJSON('gee_ex.json',flatten=TRUE)
            geedf<-as.data.frame(gee)
            colnames(geedf)<-names(gee$aoi_subregions)
            geedf$aoi_id<-c(1:nrow(geedf))
            
            
            #Intervention parameters
            j <- jsonlite::fromJSON('croplandex.json', flatten=TRUE)
            interv_sub<- j$intervention_subcategory #intervention type
            common<-as.data.frame(j$scenarios$common)
            
            bau<-as.data.frame(j$scenarios$business_as_usual)
            bau<-bau[,!names(bau)=="aoi_subregions"]
            
            bau_aoi<-as.data.frame(j$scenarios$business_as_usual$aoi_subregions)
            bau<-rbind(bau,bau)
            bau<-cbind(bau,bau_aoi)
            
            bau$scenario<-"business-as-usual"
            
            int<-as.data.frame(j$scenarios$intervention)
            int<-int[,!names(int)=="aoi_subregions"]
            int_aoi<-as.data.frame(j$scenarios$intervention$aoi_subregions)
            int<-rbind(int,int)
            int<-cbind(int,int_aoi)
            int$scenario<-"intervention"
            
            df<-rbind(bau,int)
            
            common<-rbind(common,common,common,common)
            
            df<-cbind(common,df)
            
            df <- df %>% dplyr::mutate_at(-c(13:15,29:31,34,39,40,53), as.numeric)
            
            names(df)[names(df)=="SOC_type"]<-"soil_class"
            head(df)
            #add aoi ids
            df<-merge(df,geedf[,names(geedf) %in% c("climate_zone_category","soil_class","aoi_id")],
                      by=c("climate_zone_category","soil_class"),all.x=TRUE)
      #convert uncertainty values to sd
            #Uncertainty presented as 95% CI
            df$ef_direct_n2o_sd=(df$ef_direct_n2o_uncertainty_high-df$ef_direct_n2o_uncertainty_low)/3.92
            df$ef_leaching_sd=(df$ef_leaching_uncertainty_high-df$ef_leaching_uncertainty_low)/3.92
            df$frac_org_fert_sd=(df$frac_org_fert_uncertainty_high-df$frac_org_fert_uncertainty_low)/3.92
            df$frac_leach_sd=(df$frac_leach_uncertainty_high-df$frac_leach_uncertainty_low)/3.92
            df$frac_syn_fert_sd=(df$frac_syn_fert_uncertainty_high-df$frac_syn_fert_uncertainty_low)/3.92
            #Uncertainty presented as 95% CI as percentage of mean
            df$SOC_ref_tonnes_C_ha_sd=(df$high_activity_clay_soils_HAC_error_positive/100*df$SOC_ref_tonnes_C_ha)/1.96
            #Uncertainty presented as 2 sd as percentage of mean
            df$FLU_sd=(df$FLU_error_positive/100*df$FLU)/2
            df$FMG_sd=(df$FMG_error_positive/100*df$FMG)/2
            df$FI_sd=(df$FI_error_positive/100*df$FI)/2
            #Remove other uncertainty columns 
            df<-df[,-grep("uncertainty|error",names(df))]
      #Parameters
            AOIs<-unique(df$aoi_id)
            
            
            # it be good to make paramateters like this that we might theoretically want to change 
            # more obvious, maybe bringing them to the top of the file? 
            # Just to distinquish them from the methodology driven numbers. 
            nx=500 #Number of Monte Carlo iterations
          
####GHG Calculations####
  GHGcalc<-function(i,df,nx #number of MC iterations
                      ){
        temp_bau<-df[df$aoi_id==i & df$scenario=="business-as-usual",]
        temp_int<-df[df$aoi_id==i & df$scenario=="intervention",]
        results.unc <- matrix(0, nx, 5)
        str(results.unc)
        results.unc[, 1] <- i
        results.unc[, 2] <-temp_bau$area
        results.unc[, 3] <- c(1:nx)
        colnames(results.unc) <- c('AOI','Area','Rep','SOC','N2O')
        for(m in c(1:nx)){
          ####SOC Stock Change####
              SOCREF<-rnorm(1, temp_bau$SOC_ref_tonnes_C_ha, temp_bau$SOC_ref_tonnes_C_ha_sd)
              FLUbau<-rnorm(1, temp_bau$FLU, temp_bau$FLU_sd)
              FMGbau<-rnorm(1, temp_bau$FMG, temp_bau$FMG_sd)
              FIbau<-rnorm(1, temp_bau$FI, temp_bau$FI_sd)
              FLUint<-rnorm(1, temp_int$FLU, temp_int$FLU_sd)
              FMGint<-rnorm(1, temp_int$FMG, temp_int$FMG_sd)
              FIint<-rnorm(1, temp_int$FI, temp_int$FI_sd)
              SOCbau<-SOCREF*FLUbau*FMGbau*FIbau
              SOCint<-SOCREF*FLUint*FMGint*FIint
              dSOC<-SOCint-SOCbau
              results.unc[m,4]<-dSOC*44/12
            ####N2O emissions####
                #Calculate change in N sources
                    FSN<-ifelse(grepl("nutrient management", interv_sub), #FSN - only needed for N management interventions
                                  temp_int$n_fertilizer_amount*temp_int$n_fertilizer_percent/100-
                                    temp_bau$n_fertilizer_amount*temp_bau$n_fertilizer_percent/100,0) 
                    FON<-0
                    # FON<-ifelse(grepl("nutrient management", interv_sub), ####Need separate organic N and synthetic N for volatilization if change in org amendments
                    #       temp_int$orgrate*temp_int$orgNcont/100-  ####Missing from input list - for now define as 0
                    #         temp_bau$orgrate*temp_bau$orgNcont/100, 0)
                    FSOM<-dSOC/10*1000 #FSOM 
                #Calculate N emissions
                    EFdir<-rnorm(1,temp_bau$ef_direct_n2o,temp_bau$ef_direct_n2o_sd)
                    FRAC_GASF<-.11 ####FRAC_GASF is missing from parameter list
                    FRAC_GASM<-rnorm(1,temp_bau$frac_org_fert,temp_bau$frac_org_fert_sd)
                    EF_vol<-rnorm(1,temp_bau$frac_syn_fert,temp_bau$frac_syn_fert_sd)
                    FRAC_LEACH<-rnorm(1, temp_bau$frac_leach, temp_bau$frac_leach_sd)
                    EF_leach<-rnorm(1,temp_bau$ef_leaching, temp_bau$ef_leaching_sd)
                    dirN<-(FSN+FON+FSOM)*EFdir #Direct emissions
                    volN<-(FSN*FRAC_GASF+ FON*FRAC_GASM)*EF_vol #Indirect volatilization
                    leachN<-(FSN+FON+FSOM)*FRAC_LEACH*EF_leach    #indirect leaching
                    results.unc[m,5]<-(dirN+volN+leachN)*44/28 #sum and convert to tN2O
        }
      return(results.unc)
  }
  #Set up parallel computing
      num_cores <- detectCores()
      cl <- makeCluster(num_cores)
      plan(multisession, workers = num_cores)
      clusterEvalQ(cl, library(plyr))
      clusterExport(cl, varlist = c("df"))
  #Loop through sub AOIs
      set.seed(1)
      mc<-future_lapply(
            1:length(AOIs), 
            function (x) GHGcalc(x, df, nx),
            future.seed=TRUE
          )
  #Extract mean and sd from each subAOI and store in results df 
      resdf <- matrix(0, length(AOIs), 10)
      resdf[, 1] <- AOIs
      colnames(resdf) <- c('AOI','Area','CO2_thayr','N2O_thayr','sdCO2ha','sdN2Oha',
                           'CO2_tyr','N2O_tyr','sdCO2','sdN2O')
      for(a in c(1:length(AOIs))){
        resdf[a, 2] <-mean(mc[[a]][,2])
        resdf[a, 3] <- mean(mc[[a]][,4])
        resdf[a, 4] <- mean(mc[[a]][,5])
        resdf[a, 5] <- sd(mc[[a]][,4])
        resdf[a, 6] <- sd(mc[[a]][,5])
        resdf[a, 7] <- resdf[a,2]*resdf[a,3]
        resdf[a, 8] <- resdf[a,2]*resdf[a,4]
        resdf[a, 9] <- resdf[a,2]*resdf[a,5]
        resdf[a, 10] <- resdf[a,2]*resdf[a,6]
      }
      resdf<-as.data.frame(resdf[,-2])
# convert to json output
      resdf<-resdf %>%
        group_by(AOI) %>%
        nest_legacy(.key = "Var")
      
      # this doesn't seem technically consistent with 
      # https://mathematicampr.atlassian.net/wiki/spaces/AFOLU/pages/3212444049/Code+requirements#Expected-output
      # but that might be alright if it's cool with Jenny? 
      # We need an agreed upon output format for your code though to give to the app team 
      # so we should update this or the confluence so they're consistent
      
      json_output<-jsonlite::toJSON(resdf, pretty=TRUE)
      write(json_output, paste0("output/",interv_sub,".json"))
      