# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: Grassland.R
# Author: Lisa Eash
# Date Started: 06/29/2024
# Last Edited: 06/29/2024
# Purpose: AFOLU GHG Calculations for Grassland Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: User inputs depend on subintervention category (nutrient management, change in fire
    #   management, change in livestock type or stocking rate, improved plant species, livestock grazing intensity). 
    #   -All subinterventions require:
    #       -grassland condition (nominally degraded or native, high intensity grazing, severely degraded, improved grassland)
    #       -number of improvements over initial scenario (0, 1 or 2)
    #   -For nutrient management, user must provide:
    #       -synthetic N application rate and N%. I
    #       -organic amendment application rate and N%.
    #   -For change in livestock type and stocking rate, user must provide
    #       -livestock type
    #       -stocking rate (head/ha)
    #       -average live weight of liveestock type (kg/head
    #   -For fire management, user must provide: 
    #       -Frequency of burn (number of years)
    #       -Burn timing (early or mid to late dry season)

#### Parameters and Paths ####
  # Calculations require the following parameters and their associated uncertainty:
      #SOC
          # SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
          # FLU: Land use emissions factor (always 1 for grasslands)
          # FMG: Management factor for both business as usual and intervention scenario 
          # FI: C input factor for both business as usual and intervention scenario 
      #N2O
          # FSN: Amount of synthetic N applied annually (kg N/ha/yr)
          # FON: Amount of organic N applied annually (kg N/ha/yr)
          # FPRP: Amount of N from grazing livestock waste (kg N/ha/yr)
          # EFDir: Emissions factor for direct N2O emissions from N inputs (kg N2O-N/kg N)
          # EF_PRP: Emissions factor for direct N2O emissions from grazing livestock waste (kg N2O-N/kg N)
          # FRAC_GASF: Fraction of synthetic fertilizer that volatilizes (kg NH3-N + NOx-N)/kg N
          # FRAC_GASM: Fraction of organic N that volatilizes (kg NH3-N + NOx-N)/kg N
          # EF_vol: Emissions factor for indirect N2O emissions from volatilization kg N2O-N/(kg NH3-N + NOx-N)
          # FRAC_LEACH: Fraction of N inputs that is lost to leaching (kg N loss/kg N)
          # EF_leach: Emissions factor for indirect N2O emissions from volatilization (kg N2O-N/kg N loss)
          # Nex: Nitrogen excretion rate from livestock (kg N/1000 kg animal mass/day)
          # burning_n2o_ef: Emissions factor for N2O from burning (g N2O/kg dry matter burnt)
          # combustion_factor: combustion factor for fire
          # fuel_biomass: amount of biomass available for burning (t dry matter/ha)
      #CH4
          # EF_ch4: Emissions factor for CH4 from enteric fermentation from livestock incorporation (kg CH4/head/yr)
          # burning_ch4_ef: Emissions factor for CH4 from burning (g CH4/kg dry matter burnt)

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

    # Libraries and working directory
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

    # Extract data from json input
        #Intervention parameters
            j <- jsonlite::fromJSON('grassland_ex.json', flatten=TRUE)
            interv_sub<- j$intervention_subcategory #intervention type
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
            df<-rbind(bau,int)
            common<-rbind(common,common,common,common)
            df<-cbind(common,df)
            str(df)
            df <- type.convert(df, as.is = T)
            AOIs<-unique(df$aoi_id)
        #convert uncertainty values to sd
            #Uncertainty presented as 95% CI
            unc_low<-grep("uncertainty_low",colnames(df))
            unc_low_name<-colnames(df)[unc_low]
            unc_cols<-gsub("uncertainty_low*.","\\1",unc_low_name)
            for(a in unc_cols){
              df$x<-(df[,names(df)==paste0(a,"uncertainty_high")]-
                       df[,names(df)==paste0(a,"uncertainty_low")])/3.92
              names(df)[names(df)=="x"]<-paste0(a,"sd")
            }
            df<-df[,-grep("uncertainty",names(df))]
            #Uncertainty presented as 95% CI as percentage of mean
            df$SOC_ref_tonnes_C_ha_sd=(df$high_activity_clay_soils_HAC_error_positive/100*df$SOC_ref_tonnes_C_ha)/1.96
            str(df)
            #Uncertainty presented as 2 sd as percentage of mean
            df$FMG_sd=(df$FMG_error_positive/100*df$FMG)/2
            df$FI_sd=(df$FI_error_positive/100*df$FI)/2
            #Remove other uncertainty columns 
            df<-df[,-grep("error",names(df))]


####GHG Calculations####
    GHGcalc<-function(i,df,nx #number of MC iterations
    ){
      temp_bau<-df[df$aoi_id==i & df$scenario=="business-as-usual",]
      temp_int<-df[df$aoi_id==i & df$scenario=="intervention",]
      results.unc <- matrix(0, nx, 44)
      str(results.unc)
      results.unc[, 1] <- i
      results.unc[, 2] <-temp_bau$area
      results.unc[, 3] <- c(1:nx)
      colnames(results.unc) <- c('AOI','Area','Rep','SOC',paste0('N2O_y',c(1:20)),
                                 paste0('CH4_y',c(1:20)))
      for(m in c(1:nx)){
        ####SOC Stock Change####
            SOCREF<-rnorm(1, temp_bau$SOC_ref_tonnes_C_ha, temp_bau$SOC_ref_tonnes_C_ha_sd)
            FLUbau<-temp_bau$FLU
            FMGbau<-rnorm(1, temp_bau$FMG, temp_bau$FMG_sd)
            FIbau<-rnorm(1, temp_bau$FI, temp_bau$FI_sd)
            FLUint<-temp_int$FLU
            FMGint<-rnorm(1, temp_int$FMG, temp_int$FMG_sd)
            FIint<-rnorm(1, temp_int$FI, temp_int$FI_sd)
            SOCbau<-SOCREF*FLUbau*FMGbau*FIbau
            SOCint<-SOCREF*FLUint*FMGint*FIint
            dSOC<-SOCbau-SOCint
            results.unc[m,4]<-dSOC*44/12
        ####N2O emissions####
            #Calculate change in N sources
                FSN<-ifelse(grepl("nutrient management", interv_sub),
                            temp_int$n_fertilizer_amount*temp_int$n_fertilizer_percent/100-
                              temp_bau$n_fertilizer_amount*temp_bau$n_fertilizer_percent/100,0) 
                FON<-ifelse(grepl("nutrient management", interv_sub), 
                            temp_int$org_amend_rate*temp_int$org_amend_npercent/100-  
                              temp_bau$org_amend_rate*temp_bau$org_amend_npercent/100, 0)
                FSOM<-dSOC/10*1000 #FSOM 
                if(grepl("change in livestock type or stocking rate", interv_sub)){
                  FPRP<-(temp_int$live_weight*temp_int$Nex*temp_int$stocking_rate-
                    temp_bau$live_weight*temp_bau$Nex*temp_bau$stocking_rate)/1000*365
                } else{ FPRP<-0 }
            #Calculate N emissions
                EFdir<-rnorm(1,temp_bau$ef_direct_n2o,temp_bau$ef_direct_n2o_sd)
                FRAC_GASF<-rnorm(1,temp_bau$frac_syn_fert,temp_bau$frac_syn_fert_sd)
                FRAC_GASM<-rnorm(1,temp_bau$frac_org_fert,temp_bau$frac_org_fert_sd)
                EF_vol<-rnorm(1,temp_bau$ef_vol,temp_bau$ef_vol_sd)
                FRAC_LEACH<-rnorm(1, temp_bau$frac_leach, temp_bau$frac_leach_sd)
                EF_leach<-rnorm(1,temp_bau$ef_leaching, temp_bau$ef_leaching_sd)
                EF_PRP<-rnorm(1, temp_bau$ef_prp, temp_bau$ef_prp_sd)
                dirN<-(FSN+FON+FSOM)*EFdir + FPRP*EF_PRP #Direct emissions
                volN<-(FSN*FRAC_GASF+ (FON+FPRP)*FRAC_GASM)*EF_vol #Indirect volatilization
                leachN<-(FSN+FON+FSOM+FPRP)*FRAC_LEACH*EF_leach    #indirect leaching
                N2O<-(dirN+volN+leachN)/1000*44/28 #sum and convert to tN2O
                results.unc[m,c(5:24)]<-N2O
                #Add change in N2O due in fire management
                if(grepl("change in fire management", interv_sub)){
                  fire_n2o_ef<-rnorm(1,temp_bau$burning_n2o_ef_mean,temp_bau$burning_n2o_ef_sd)
                  #Add fires for bau scenario 
                  if(temp_bau$fire_used=="True"){
                    CF_bau<-rnorm(1,temp_bau$combustion_factor_mean,temp_bau$combustion_factor_sd)
                    MB_bau<-rnorm(1,temp_bau$fuel_biomass_mean,temp_bau$fuel_biomass_sd)
                    fireN2O_bau<-MB_bau*CF_bau*fire_n2o_ef/1000
                    results.unc[m,5]<-N2O-fireN2O_bau
                    fire_per_bau<-temp_bau$fire_management_years
                    fire_yrs_bau<-vector()
                    for (y in 1:19) { 
                      if (y%%fire_per_bau == 0) {
                        fire_yrs_bau<-append(fire_yrs_bau,y)
                      }
                    }
                    results.unc[m,5+fire_yrs_bau]<-N2O-fireN2O_bau
                  }
                  if(temp_int$fire_used=="True"){
                    CF_int<-rnorm(1,temp_int$combustion_factor_mean,temp_int$combustion_factor_sd)
                    MB_int<-rnorm(1,temp_int$fuel_biomass_mean,temp_int$fuel_biomass_sd)
                    fireN2O_int<-MB_int*CF_int*fire_n2o_ef/1000
                    results.unc[m,5]<-results.unc[m,5]+fireN2O_int
                    fire_per_int<-temp_int$fire_management_years
                    fire_yrs_int<-vector()
                    for (y in 1:19) { 
                      if (y%%fire_per_int == 0) {
                        fire_yrs_int<-append(fire_yrs_int,y)
                      }
                    }
                    results.unc[m,5+fire_yrs_int]<-results.unc[m,5+fire_yrs_int]+fireN2O_int
                  }
                } 
                
        ####CH4 emissions####
            #emissions from enteric fermentation (if livestock incorporation)
            if(grepl("change in livestock type or stocking rate", interv_sub)){
              EF_CH4<-rnorm(1, temp_int$ef_ch4, 
                            temp_int$ef_ch4*.50/2) #In Section 10.3.4 IPCC (2006), uncertainty range (2 sd) is 30-50% for livestock parameters
              CH4_live<-(temp_int$stocking_rate-temp_bau$stocking_rate)*EF_CH4/1000
            } else{ CH4_live<-0 }
            results.unc[m,c(25:44)]<-CH4_live
            #emissions from fire management
            if(grepl("change in fire management", interv_sub)){
              fire_ch4_ef<-rnorm(1,temp_bau$burning_ch4_ef_mean,temp_bau$burning_ch4_ef_sd)
              #Add fires for bau scenario 
              if(temp_bau$fire_used=="True"){
                fireCH4_bau<-MB_bau*CF_bau*fire_ch4_ef/1000
                results.unc[m,25]<-CH4_live-fireCH4_bau
                results.unc[m,25+fire_yrs_bau]<-CH4_live-fireCH4_bau
              }
              if(temp_int$fire_used=="True"){
                fireCH4_int<-MB_int*CF_int*fire_ch4_ef/1000
                results.unc[m,25]<-results.unc[m,25]+fireCH4_int
                results.unc[m,25+fire_yrs_int]<-results.unc[m,25+fire_yrs_int]+fireCH4_int
              }
              } 
      }
      results.unc<-lapply(as.data.frame(results.unc), as.numeric)
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
        resdf <- matrix(0, length(AOIs), 166)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area','CO2_thayr',
                             paste0('N2O_thayr_y',c(1:20)),
                             paste0('CH4_thayr_y',c(1:20)),
                             'sdCO2ha',paste0('sdN2Oha_y',c(1:20)),
                             paste0('sdCH4ha_y',c(1:20)),
                             'CO2_tyr',paste0('N2O_tyr_y',c(1:20)),
                             paste0('CH4_tyr_y',c(1:20)),
                             'sdCO2',
                             paste0('sdN2O_y',c(1:20)),
                             paste0('sdCH4_y',c(1:20)))
        for(a in c(1:length(AOIs))){
          temp_res<-as.data.frame(mc[[a]])
          resdf[a, 2] <-mean(temp_res$Area)
          for(z in c(4:44)){
            resdf[a, z-1] <- mean(temp_res[,z])
            resdf[a, z+40] <- sd(temp_res[,z])
          }
        }
        for(k in c(3:84)){
          resdf[, 82+k] <- resdf[,2]*resdf[,k]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        str(resdf)
        restib<-data.frame()
        for(r in c(1:nrow(resdf))){
          aoi_tib <-tibble(
            AOI = paste0("Sub_AOI-",r),
            CO2_thayr = list(rep(resdf[r,2],20)),
            sdCO2ha = list(rep(resdf[r,43],20)),
            CO2_tyr = list(rep(resdf[r,84],20)),
            sdCO2 = list(rep(resdf[r,125],20)),
            N2O_thayr = list(unname(unlist(resdf[r,c(3:22)]))),
            sdN2Oha = list(unname(unlist(resdf[r,c(44:63)]))),
            N2O_tyr = list(unname(unlist(resdf[r,c(85:104)]))),
            sdN2O = list(unname(unlist(resdf[r,c(126:145)]))),
            CH4_thayr = list(unname(unlist(resdf[r,c(23:42)]))),
            sdCH4ha  = list(unname(unlist(resdf[r,c(64:83)]))),
            CH4_tyr = list(unname(unlist(resdf[r,c(105:124)]))),
            sdCH4  = list(unname(unlist(resdf[r,c(146:165)])))
          )
          restib<-rbind(restib,aoi_tib)
        }
        json_output<-jsonlite::toJSON(restib,
                                      auto_unbox = T,
                                      pretty = T)
        write(json_output, paste0("output/",interv_sub,".json"))
