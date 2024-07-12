# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: Agroforestry.R
# Author: Lisa Eash
# Date Started: 07/09/2024
# Last Edited: 07/11/2024
# Purpose: AFOLU GHG Calculations for Agroforestry Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: The user must select the subintervention of 'adding woody biomass to cropland' or 
    #   'adding woody biomass to grassland'. 
    #   -Inputs for all subinterventions for both the business-as-usual and intervention scenarios are: 
    #       -synthetic N application rate (kg/ha) and N%
    #       -organic amendment application rate (kg/ha) and N%
    #       -type of agroforestry (only for intervention scenario)
    #       -portion of area converted to agroforestry (%)
    #   -Inputs for 'adding woody biomass to grassland' are:
    #       -grassland condition (nominally degraded or native, high intensity grazing, severely degraded, improved grassland)
    #       -number of improvements over initial scenario (0, 1 or 2)
    #   -Inputs for the cropland scenario are: 
    #       -crop type (annual or perennial)
    #       -organic matter management (low, medium, high, high with manure)
    #       -tillage type (full, reduced, no-till, or NA)

#### Parameters and Paths ####
  # Calculation requires the following parameters and their associated uncertainty:
      #Biomass C
          # crop_biomass: Amount of crop biomass lost or gained 1 year after conversion (t C/ha; only if adding woody biomass to cropland)
          # grass_biomass: Amount of grassland biomass lost or gained 1 year after conversion (t dry matter/ha; only if adding woody biomass to grassland)
          # biomass_c_rate: Annual C growth rate of agroforestry system (t C/ha/year)
          # max_biomass_c: Max C biomass of agroforestry system (t C/ha)
      #SOC
          # SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
          # FLU: Land use emissions factor (always 1 for grasslands)
          # FMG: Management factor for both business as usual and intervention scenario 
          # FI: C input factor for both business as usual and intervention scenario 
      #N2O
          # EFDir: Emissions factor for direct N2O emissions from N inputs (kg N2O-N/kg N)
          # FRAC_GASF: Fraction of synthetic fertilizer that volatilizes (kg NH3-N + NOx-N)/kg N
          # FRAC_GASM: Fraction of organic N that volatilizes (kg NH3-N + NOx-N)/kg N
          # EF_vol: Emissions factor for indirect N2O emissions from volatilization kg N2O-N/(kg NH3-N + NOx-N)
          # FRAC_LEACH: Fraction of N inputs that is lost to leaching (kg N loss/kg N)
          # EF_leach: Emissions factor for indirect N2O emissions from volatilization (kg N2O-N/kg N loss)
      #CH4- assumed to be 0 

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
    scen=c("business_as_usual","intervention")
    Cf_herbaceous <- 0.47 #IPCC (2006)
    FLU_forest<-1
    FI_forest<-1
    FMG_forest<-1

    # Extract data from json input
        #Intervention parameters
            j <- jsonlite::fromJSON('agroforestry_ex.json', flatten=TRUE)
            interv_sub<- j$intervention_subcategory #intervention type
            nyears<-as.numeric(j$intervention_years_to_predict)
            common<-as.data.frame(j$scenarios$common)
            df_bau<-as.data.frame(j$scenarios$business_as_usual)
            df_bau<-df_bau[,!names(df_bau)=="aoi_subregions"]
            df_bau_aoi<-as.data.frame(j$scenarios$business_as_usual$aoi_subregions)
            df_bau_aoi$aoi_id<-c(1:nrow(df_bau_aoi))
            df_bau<-rbind(df_bau,df_bau)
            df_bau<-cbind(df_bau,df_bau_aoi)
            df_bau$scenario<-"business_as_usual"
            df_int<-as.data.frame(j$scenarios$intervention)
            df_int<-df_int[,!names(df_int)=="aoi_subregions"]
            df_int_aoi<-as.data.frame(j$scenarios$intervention$aoi_subregions)
            df_int_aoi$aoi_id<-c(1:nrow(df_int_aoi))
            df_int<-rbind(df_int,df_int)
            df_int<-cbind(df_int,df_int_aoi)
            df_int$scenario<-"intervention"
            common<-rbind(common,common)
            df_int<-cbind(common,df_int)
            df_bau<-cbind(common,df_bau)
            df_int <- type.convert(df_int, as.is = T)
            df_bau <- type.convert(df_bau, as.is = T)
            AOIs<-unique(df_int$aoi_id)
        #convert uncertainty values to sd
          for(s in scen){
            if(s=="business_as_usual"){
              df=df_bau
            } else {df=df_int}
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
            df<-df[,-grep("high_activity_clay_soils_HAC_error_",names(df))]
            #Uncertainty presented as 2 sd as percentage of mean
            err_low<-grep("error_negative",colnames(df))
            err_low_name<-colnames(df)[err_low]
            err_cols<-gsub("_error_negative*.","\\1",err_low_name)
            for(a in err_cols){
              df$x<-(df[,names(df)==paste0(a,"_error_positive")]/100*
                       df[,names(df)==a])/2
              names(df)[names(df)=="x"]<-paste0(a,"_sd")
            }
            #Remove other uncertainty columns 
            df<-df[,-grep("error",names(df))]
            if(s=="business_as_usual"){
              df_bau=df
            } else {df_int=df}
          }

####GHG Calculations####
    GHGcalc<-function(i,df,nx #number of MC iterations
    ){
      temp_bau<-df_bau[df_bau$aoi_id==i,]
      temp_int<-df_int[df_int$aoi_id==i,]
      agroforestry_type<-temp_int$agroforestry_type
      results.unc <- matrix(0, nx, 3+nyears+2)
      str(results.unc)
      results.unc[, 1] <- i
      results.unc[, 2] <-temp_bau$area
      results.unc[, 3] <- c(1:nx)
      colnames(results.unc) <- c('AOI','Area','Rep',paste0('CO2_y',c(1:20)),'N2O',
                                 'CH4')
      for(m in c(1:nx)){
        ####SOC Stock Change####
                SOCREF<-rnorm(1, temp_bau$SOC_ref_tonnes_C_ha, temp_bau$SOC_ref_tonnes_C_ha_sd) #same for grass and crop
            #Calculate SOC of bau scenario
                FLU<-rnorm(1, temp_bau$FLU, temp_bau$FLU_sd)
                FMG_bau<-rnorm(1, temp_bau$FMG, temp_bau$FMG_sd)
                FI_bau<-rnorm(1, temp_bau$FI, temp_bau$FI_sd)
                SOC_bau<-SOCREF*FLU*FMG_bau*FI_bau
            #Calculate SOC of intervention on non-forested land
                if(temp_int$agroforestry_land_portion!=100){
                FLU_int<-rnorm(1, temp_bau$FLU, temp_bau$FLU_sd)
                FMG_int<-rnorm(1, temp_bau$FMG, temp_bau$FMG_sd)
                FI_int<-rnorm(1, temp_bau$FI, temp_bau$FI_sd)
                SOC_nonforest<-SOCREF*FLU*FMG_int*FI_int
                } else {
                  SOC_nonforest<-SOC_bau
                }
                SOC_forest<-SOCREF*FLU_forest*FMG_forest*FI_forest
                SOC_int<-(SOC_forest*temp_int$agroforestry_land_portion/100+
                            SOC_nonforest*(1-(temp_int$agroforestry_land_portion/100)))
            #Calculate difference
                  dSOC<-(SOC_bau-SOC_int)/20*44/12
        ####Biomass C#### 
            #Calculate gain in annual agroforestry biomass 
                BioTree<-rnorm(1,temp_int$biomass_c_rate,temp_int$biomass_c_rate_sd)*44/12*temp_int$agroforestry_land_portion/100
                MaxBioTree<-rnorm(1,temp_int$max_biomass_c,temp_int$max_biomass_c_sd)*44/12*temp_int$agroforestry_land_portion/100
                max_year<-ceiling(MaxBioTree/BioTree)
                BioTree_max_year<-MaxBioTree-(max_year-1)*BioTree
            #Calculate loss of initial biomass in grassland or cropland
                if(interv_sub=="Adding woody crops to cropland"){
                  bioCloss<-rnorm(1, temp_bau$crop_biomass, temp_bau$crop_biomass_sd)*44/12
                } else{
                  grass_bio<-rnorm(1,temp_bau$grass_biomass,temp_bau$grass_biomass_sd)
                  bioCloss<-Cf_herbaceous*grass_bio*44/12*temp_int$agroforestry_land_portion/100
                }
            #Store yearly change in biomass and SOC change in results matrix
                results.unc[m,c(3+1)]<-dSOC-BioTree+bioCloss #year 1
                results.unc[m,c(3+c(2:min(20,(max_year-1))))]<-dSOC-BioTree #prior to max agroforestry biomass
                if(max_year<=20){
                  results.unc[m,c(3+max_year)]<-dSOC-BioTree_max_year #when agroforestry reaches max biomass
                  if(max_year<20){
                  results.unc[m,c(3+c((max_year+1):nyears))]<-dSOC #after agroforestry reaches max biomass
                  }
                }
        ####N2O emissions####
            #Calculate change in N sources- N inputs in grassland vs. N inputs in cropland, assumed to be constant
                FSN<-(temp_int$n_fertilizer_amount*temp_int$n_fertilizer_percent/100)-
                  (temp_bau$n_fertilizer_amount*temp_bau$n_fertilizer_percent/100)
                FON<-(temp_int$org_amend_rate*temp_int$org_amend_npercent/100)-
                  (temp_bau$org_amend_rate*temp_bau$org_amend_npercent/100)
                FSOM<-dSOC/10*1000 #FSOM 
            #Calculate N emissions
                EFdir<-rnorm(1,temp_bau$ef_direct_n2o,temp_bau$ef_direct_n2o_sd)
                FRAC_GASF<-rnorm(1,temp_bau$frac_syn_fert,temp_bau$frac_syn_fert_sd)
                FRAC_GASM<-rnorm(1,temp_bau$frac_org_fert,temp_bau$frac_org_fert_sd)
                EF_vol<-rnorm(1,temp_bau$ef_vol,temp_bau$ef_vol_sd)
                FRAC_LEACH<-rnorm(1, temp_bau$frac_leach, temp_bau$frac_leach_sd)
                EF_leach<-rnorm(1,temp_bau$ef_leaching, temp_bau$ef_leaching_sd)
                dirN<-(FSN+FON+FSOM)*EFdir #Direct emissions
                volN<-(FSN*FRAC_GASF+ FON*FRAC_GASM)*EF_vol #Indirect volatilization
                leachN<-(FSN+FON+FSOM)*FRAC_LEACH*EF_leach    #indirect leaching
                N2O<-(dirN+volN+leachN)/1000*44/28 #sum and convert to tN2O
                results.unc[m,c(24)]<-N2O
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
        resdf <- matrix(0, length(AOIs), 2+nyears*4+2*4)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area',paste0('CO2_thayr_y',c(1:20)),
                             'N2O_thayr','CH4_thayr',
                             paste0('sdCO2ha_y',c(1:20)),
                             'sdN2Oha','sdCH4ha',
                             paste0('CO2_tyr_y',c(1:20)),
                             'N2O_tyr','CH4_tyr',
                             paste0('sdCO2_y',c(1:20)),
                             'sdN2O','sdCH4')
        for(a in c(1:length(AOIs))){
          temp_res<-as.data.frame(mc[[a]])
          resdf[a, 2] <-mean(temp_res$Area)
          for(z in c(4:25)){
            resdf[a, z-1] <- mean(temp_res[,z])
            resdf[a, z+21] <- sd(temp_res[,z])
          }
        }
        for(k in c(3:46)){
          resdf[, 44+k] <- resdf[,2]*resdf[,k]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        restib<-data.frame()
        for(r in c(1:nrow(resdf))){
          aoi_tib <-tibble(
            AOI = paste0("Sub_AOI-",r),
            CO2_thayr = list(unname(unlist(resdf[r,c(2:21)]))),
            sdCO2ha = list(unname(unlist(resdf[r,c(24:43)]))),
            CO2_tyr = list(unname(unlist(resdf[r,c(46:65)]))),
            sdCO2 = list(unname(unlist(resdf[r,c(68:87)]))),
            N2O_thayr = list(rep(resdf[r,22],20)),
            sdN2Oha = list(rep(resdf[r,44],20)),
            N2O_tyr = list(rep(resdf[r,66],20)),
            sdN2O = list(rep(resdf[r,88],20)),
            CH4_thayr = list(rep(resdf[r,23],20)),
            sdCH4ha  = list(rep(resdf[r,45],20)),
            CH4_tyr = list(rep(resdf[r,67],20)),
            sdCH4  = list(rep(resdf[r,89],20))
          )
          restib<-rbind(restib,aoi_tib)
        }
        json_output<-jsonlite::toJSON(restib,
                                      auto_unbox = T,
                                      pretty = T)
        write(json_output, paste0("output/",interv_sub,".json"))
