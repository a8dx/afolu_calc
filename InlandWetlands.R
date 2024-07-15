# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: InlandWetlands.R
# Author: Lisa Eash
# Date Started: 07/01/2024
# Last Edited: 07/03/2024
# Purpose: AFOLU GHG Calculations for Inland Wetlands Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: The user must first select all intervention subtypes that apply (draining,
    #   rewetting, peat extraction, land use conversion). Other user inputs include the following:
    #       -wetland_area: Area of wetland included in project (ha)
    #       -land_use_type_init: Initial land use type prior to project start (Natural Forest, Forest Plantation, Cropland, Rice, Grassland, Peatland, Settlements, Other land)
    #       -plantation_type: If Forest Plantation selected for land_use_type_init, type of plantation
    #       -fire_used: Whether or not fire was used for land conversion
    #       -fire_frequency: expected frequency of peatland fires throughout project
    #       -peat_extraction: If intervention subtype includes peat extraction, annual weight of peat extraction (tonnes/year)
    #       -nutrient_status: Wetland nutrient status at start of project (poor, rich)


#### Parameters and Paths ####
  # Calculations requires the following parameters and their associated uncertainty:
      #Emissions from soil draining
          # EF_CO2_drained: CO2 emissions factor from soil draining (t CO2-C/ha/yr)
          # doc_flux_natural: natural DOC flux for organic soils (t C/ha/yr)
          # dDOC: proportional increase in DOC due to drainage
          # frac_doc: conversion factor for proportion of DOC converted to CO2 following export from site 
          # EF_CH4: CH4 emissions factor from soil draining (kg CH4/ha/yr)
          # EF_N2O: N2O emissions factor from soil draining (kg N2O-N/ha/yr)
      #Emissions from soil rewetting
          # EF_CO2_rewet: CO2 emissions factor from soil rewetting (t CO2-C/ha/yr)
          # EF_DOC_rewet: DOC emissions factor from soil rewetting (t CO2-C/ha/yr)
          # EF_CH4_rewet: CH4 emissions factor from soil rewetting (kg CH4-C/ha/yr)
      #Emissions from peat extraction
          # cfraction_peat: carbon fraction of peat extracted (t C/t peat)
          # EF_CO2_peat: CO2 emissions factor for lands managed for peat extraction (t C/ha/yr)
          # EF_N2O_peat: N2O emissions factor for lands managed for peat extraction (kg N2O-N/ha/yr)
      #Emissions from fire during land use conversion
          # crop_biomass: Biomass for cropland (t C/ha; if initial land use type = cropland)
          # grass_biomass: Dead organic matter plus live biomass for grassland (t dry matter/ha; if initial land use type = grassland)
          # forest_biomass: Aboveground forest biomass for natural or plantation forests (if initial land use type; t dry matter/ha) ##Lisa, we're using variables called gedi_l4b_agb_mean_total_mg, and gedi_l4b_agb_se_total_mg, for the forest aboveground biomass and standard error. This are going to come from remote sensing and will be part of Step 1, so they'll be printed for each sub-aoi. Units are Mg C for the entire aoi, not per ha. Could you please modify the code to reflect this?##
          # forest_R: Ratio of belowground:aboveground forest biomass for natural or plantation forests (if initial land use type)
          # forest_litter_c: Litter forest C for natural or plantation forests (if initial land use type; t C/ha)
          # forest_deadom_c: Dead forest C for natural or plantation forests (if initial land use type; t C/ha)
          # combustion_factor_init: Combustion factor for fire during land use conversion given initial land use type (only if fire used in land conversion)
          # burning_co2_ef_init: CO2 emissions factor for fire during land use conversion (g CO2/kg d.m burned; only if fire used in land conversion from forests)
          # burning_n2o_ef_init: N2O emissions factor for fire during land use conversion given initial land use type (g N2O/kg d.m burned; only if fire used in land conversion)
          # burning_ch4_ef_init: CH4 emissions factor for fire during land use conversion given initial land use type (g CH4/kg d.m burned; only if fire used in land conversion)
      #Emissions from fires throughout project 
          # MB_peat: biomass available for burning on peatlands (t dry matter/ha)
          # EF_CH4_fire_peat: CH4 emissions factor for peatland fire (g CH4/kg dry matter)
          # EF_CO2_fire_peat: CO2 emissions factor for peatland fire (g CO2-C/kg dry matter)

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
    Cf_herbaceous <- 0.47 #Carbon fraction herbaceous biomass IPCC (2006)
    Cf_forest=0.47 #Carbon fraction forest biomass IPCC (2019)
    combustion_factor_peat=1 #Table 2.7
    scen=c("business_as_usual","intervention")

    # Extract data from json input
        #Intervention parameters
            j <- jsonlite::fromJSON('inlandwetlands_ex.json', flatten=TRUE)
            intervention_type<-j$intervention_type
            interv_luc<- j$intervention_luc
            nyears<-as.numeric(j$intervention_years_to_predict)
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
            if(df[1,"land_use_type"] %in% c("Natural Forest","Forest Plantation")){
              df$forest_R_sd<-NA
              for(r in c(1:nrow(df))){
                if(df[r,"forest_R_uncertainty_type"]=="default"){
                  df[r,"forest_R_sd"]<-(df[r,"forest_R_uncertainty"]/100*
                                          df[r,"forest_R"])/1.96
                }
                if(df[r,"forest_R_uncertainty_type"]=="SD"){
                  df[r,"forest_R_sd"]<-df[r,"forest_R_uncertainty"]
                }
              }
            }
            df<-df[,-grep("uncertainty",names(df))]
            #Uncertainty presented as +/- 95% CI
            unc_ci<-grep("CI_95_positive",colnames(df))
            unc_ci_name<-colnames(df)[unc_ci]
            unc_cols<-gsub("CI_95_positive*.","\\1",unc_ci_name)
            for(a in unc_cols){
              df$x<-(df[,names(df)==paste0(a,"CI_95_positive")])/1.96
              names(df)[names(df)=="x"]<-paste0(a,"sd")
            }
            #Uncertainty presented as 2 sd as percentage of mean
            err_low<-grep("error_negative",colnames(df))
            err_low_name<-colnames(df)[err_low]
            err_cols<-gsub("_error_negative*.","\\1",err_low_name)
            for(a in err_cols){
              df$x<-(df[,names(df)==paste0(a,"_error_positive")]/100*
                       df[,names(df)==a])/2
              names(df)[names(df)=="x"]<-paste0(a,"_sd")
            }
            #If emissions factor sds are NA, assume 40% uncertainty
            unc_sd<-grep("sd",colnames(df))
            unc_sd_name<-colnames(df)[unc_sd]
            for(z in unc_sd_name){
              for(r in c(1:nrow(df))){
                if(is.na(df[r,names(df)==z])){
                  df[r,names(df)==z]<-df[r,names(df)==gsub("_sd*.","\\1",z)]*.40/1.96
                }
              }
            }
            if(s=="business_as_usual"){
              df_bau=df
            } else {df_int=df}
            }


####GHG Calculations####
    GHGcalc<-function(i,df_bau,df_int,nx #number of MC iterations
    ){
      
      results.unc <- matrix(0, nx, 3+nyears*4+2*2)
      results.unc[, 1] <- i
      results.unc[, 2] <-df_bau[df_bau$aoi_id==i,]$area
      results.unc[, 3] <- c(1:nx)
      colnames(results.unc) <- c('AOI','Area','Rep',paste0('CO2_y',c(1:20),'_bau'),
                                 paste0('N2O_y',c(1:2),'_bau'),
                                 paste0('CH4_y',c(1:20),'_bau'),
                                 paste0('CO2_y',c(1:20),'_int'),
                                 paste0('N2O_y',c(1:2),'_int'),
                                 paste0('CH4_y',c(1:20),'_int'))
    for(j in c(1:length(scen))){
      s<-scen[j]
      if(s=="business_as_usual"){
        temp=df_bau[df_bau$aoi_id==i,]
      } else {temp=df_int[df_int$aoi_id==i,]}
      land_use_type<-temp$land_use_type
      for(m in c(1:nx)){
        ####Land conversion- forest biomass removal and fires####
            if(interv_luc=="yes" && s=="business_as_usual"){
              if(land_use_type %in% c("Natural Forest","Plantation Forest")){
                AGBiomass<-rnorm(1,temp$forest_biomass, temp$forest_biomass_sd)
                R<-rnorm(1,temp$forest_R, temp$forest_R_sd)
                DOM<-rnorm(1,temp$forest_deadom_c,temp$forest_deadom_c_sd)
                Litter<-rnorm(1, temp$forest_litter_c,temp$forest_litter_c_sd)
                BioC<-(1+R)*AGBiomass*Cf_forest+DOM+Litter
              } else {BioC<-0}
              if(temp$fire_used=="True"){
                CF<-rnorm(1,temp$combustion_factor_init,temp$combustion_factor_init_sd)
                fire_n2o_ef<-rnorm(1,temp$burning_n2o_ef_init,temp$burning_n2o_ef_init_sd)
                fire_ch4_ef<-rnorm(1,temp$burning_ch4_ef_init,temp$burning_ch4_ef_init_sd)
                if(land_use_type %in% c("Natural Forest","Plantation Forest")){
                  MB<-BioC/Cf_forest
                  fire_co2_ef<-rnorm(1,temp$burning_co2_ef_init,temp$burning_co2_ef_init_sd)
                  fireCO2<--MB*CF*fire_co2_ef/1000
                }else{fireCO2=0}
                if(land_use_type %in% c("Cropland","Rice")){
                  BioCrop<-rnorm(1,temp$crop_biomass,temp$crop_biomass_sd)
                  MB<-BioCrop/Cf_herbaceous
                }
                if(land_use_type == "Grassland"){
                  MB<-rnorm(1,temp$grass_biomass,temp$grass_biomass_sd)
                }
                if(land_use_type == "Peatland Extraction"){
                  MB<-rnorm(1,temp$MB_peat,temp$MB_peat_sd)
                }
                fireN2O<--MB*CF*fire_n2o_ef/1000
                fireCH4<--MB*CF*fire_ch4_ef/1000
                BioC<-0
              } else{
                BioC=BioC*44/12
                fireCO2<-0
                fireN2O<-0
                fireCH4<-0
              } 
            } else {
              BioC=0
              fireCO2<-0
              fireN2O<-0
              fireCH4<-0
            } 
            
        ####Soil emissions####
          if(land_use_type!="Rewetting for wetland restoration"){
            drain<-temp$wetland_draining
            #CO2
                C_onsite<-rnorm(1, temp$EF_CO2_drained, temp$EF_CO2_drained_sd) 
                DOC_flux_natural<-rnorm(1, temp$doc_flux_natural, temp$doc_flux_natural_sd) 
                if(drain=="yes"){ #if not drained, only account for natural doc flux. If drained, include dDOC due to draining
                  dDOC<-rnorm(1, temp$dDOC, temp$dDOC_sd)
                } else{dDOC=0}
                Frac_doc<-rnorm(1, temp$frac_doc, temp$frac_doc_sd)
                EFDOC<-DOC_flux_natural*(1+dDOC)*Frac_doc 
                CO2_soil<-(C_onsite+EFDOC)*44/12
            if(drain=="yes"){
            #CH4
                CH4_soil<-rnorm(1, temp$EF_CH4,temp$EF_CH4_sd)/1000
            #N2O
                N2O_soil<-rnorm(1, temp$EF_N2O,temp$EF_N2O_sd)*44/28/1000
            } else{
              CH4_soil<-0
              N2O_soil<-0
            }
          } else{
            #Emissions from soil rewetting (for wetlands restoration)
                #CO2
                CO2_onsite=rnorm(1, temp$EF_CO2_rewet,temp$EF_CO2_rewet_sd)
                CO2_offsite=rnorm(1, temp$EF_DOC_rewet,temp$EF_DOC_rewet_sd)
                CO2_soil<-(CO2_onsite+CO2_offsite)*44/12
                #CH4
                CH4_soil<-rnorm(1, temp$EF_CH4_rewet,temp$EF_CH4_rewet_sd)/1000
                #N2O
                N2O_soil<-0
          }
        ####Peat extraction####
              if(land_use_type=="Peat Extraction"){
              #CO2
                CO2_onsite=rnorm(1, temp$EF_CO2_peat,temp$EF_CO2_peat_sd) 
                CO2_offsite=temp$peat_extraction*temp$cfraction_peat 
                CO2_peat<-(CO2_onsite+CO2_offsite)*44/12
              #N2O
                N2O_peat<-rnorm(1, temp$EF_N2O_peat,temp$EF_N2O_peat_sd)*44/28/1000 
              } else{
                CO2_peat<-0
                N2O_peat<-0
              }
        ####Peat fires####
              if(land_use_type %in% c("Peatland Extraction","Rewetting for wetland restoration") && 
                 temp$fire_frequency>0){
                EF_CO2_fire_peat=rnorm(1, temp$EF_CO2_fire_peat, temp$EF_CO2_fire_peat_sd)
                EF_CH4_fire_peat=rnorm(1, temp$EF_CH4_fire_peat, temp$EF_CH4_fire_peat_sd)
                MB_peat=rnorm(1, temp$MB_peat, temp$MB_peat_sd)
                CO2_fire_peat<-MB_peat*combustion_factor_peat*EF_CO2_fire_peat*44/12/1000
                CH4_fire_peat<-MB_peat*combustion_factor_peat*EF_CH4_fire_peat/1000
              } else{
                CO2_fire_peat<-0
                CH4_fire_peat<-0
              }
              
              
              
    ####Record in matrix####
        results.unc[m,c(4+
                          (j-1)*(nyears*2+2))]<-CO2_soil+CO2_peat-BioC-fireCO2
        results.unc[m,c((5+(j-1)*(nyears*2+2)):
                          ((3+nyears)+(j-1)*(nyears*2+2)))]<-CO2_soil+CO2_peat
        results.unc[m,c(4+nyears+(j-1)*(nyears*2+2))]<-N2O_soil+N2O_peat-fireN2O
        results.unc[m,c(5+nyears+(j-1)*(nyears*2+2))]<-N2O_soil+N2O_peat
        results.unc[m,c(6+nyears+(j-1)*(nyears*2+2))]<-CH4_soil-fireCH4
        results.unc[m,c(((7+nyears)+(j-1)*(nyears*2+2)):
                          ((5+(nyears*2))+(j-1)*(nyears*2+2)))]<-CH4_soil
        if(land_use_type %in% c("Peatland Extraction","Rewetting for wetland restoration")){
            if(temp$fire_frequency>0){
              fire_per<-temp$fire_frequency
              fire_yrs<-vector()
              for (y in c(1:(nyears-1))) { 
                if (y%%fire_per == 0) {
                  fire_yrs<-append(fire_yrs,y+1)
                }
              }
              results.unc[m,3+fire_yrs+(j-1)*(nyears*2+2)]<-results.unc[m,3+fire_yrs+(j-1)*(nyears*2+2)]+CO2_fire_peat
              results.unc[m,5+nyears+fire_yrs+(j-1)*(nyears*2+2)]<-results.unc[m,5+nyears+fire_yrs+(j-1)*(nyears*2+2)]+CH4_fire_peat
            }
        }
      }
    }
      results.unc<-as.data.frame(results.unc)
      return(results.unc)
    }
      
    #Set up parallel computing
        num_cores <- detectCores()
        cl <- makeCluster(num_cores)
        plan(multisession, workers = num_cores)
        clusterExport(cl, varlist = c("df"))
    #Loop through sub AOIs
        mc<-future_lapply(
          1:length(AOIs), 
          function (x) GHGcalc(x, df_bau,df_int, nx),
          future.seed=TRUE
        )
    #Extract mean and sd from each subAOI and store in results df 
        resdf <- matrix(0, length(AOIs), 2+nyears*8+4*2)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area',paste0('CO2_thayr_y',c(1:20)),
                             paste0('N2O_thayr_y',c(1:2)),
                             paste0('CH4_thayr_y',c(1:20)),
                             paste0('sdCO2ha_y',c(1:20)),
                             paste0('sdN2Oha_y',c(1:2)),
                             paste0('sdCH4ha_y',c(1:20)),
                             paste0('CO2_tyr_y',c(1:20)),
                             paste0('N2O_tyr_y',c(1:2)),
                             paste0('CH4_tyr_y',c(1:20)),
                             paste0('sdCO2_y',c(1:20)),
                             paste0('sdN2O_y',c(1:2)),
                             paste0('sdCH4_y',c(1:20)))
    #Calculate difference in scenarios for all MC reps
        for(a in c(1:length(AOIs))){
          temp_res<-as.data.frame(mc[[a]])
          temp_res[,c(88:(87+(nyears*2+2)))]<-NA
          for(y in c(1:(nyears*2+2))){
            temp_res[, 87+y] <- temp_res[, 45+y]-temp_res[, 3+y]
          }
          temp_res<-temp_res[,-c(4:87)]
          colnames(temp_res)<-c('AOI','Area','Rep',paste0('CO2_y',c(1:20)),
                                'N2O_y1','N2O_y2',
                                paste0('CH4_y',c(1:20)))
          resdf[a, 2] <-mean(temp_res$Area)
          for(z in c(1:(nyears*2+2))){
            resdf[a, z+2] <- mean(temp_res[,z+3])
            resdf[a, z+2+(nyears*2+2)] <- sd(temp_res[,z+3])
          }
        }
        for(k in c(1:(nyears*4+2*2))){
          resdf[, 2+(nyears*4+2*2)+k] <- resdf[,2]*resdf[,2+k]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        restib<-data.frame()
        for(r in c(1:nrow(resdf))){
          aoi_tib <-tibble(
            AOI = paste0("Sub_AOI-",r),
            CO2_thayr = list(unname(unlist(resdf[r,1+c(1:nyears)]))),
            sdCO2ha = list(unname(unlist(resdf[r,43+c(1:nyears)]))),
            CO2_tyr = list(unname(unlist(resdf[r,85+c(1:nyears)]))),
            sdCO2 = list(unname(unlist(resdf[r,127+c(1:nyears)]))),
            N2O_thayr = list(c(resdf[r,22],rep(resdf[r,23],19))),
            sdN2Oha = list(c(resdf[r,64],rep(resdf[r,65],19))),
            N2O_tyr = list(c(resdf[r,106],rep(resdf[r,107],19))),
            sdN2O = list(c(resdf[r,148],rep(resdf[r,149],19))),
            CH4_thayr = list(unname(unlist(resdf[r,23+c(1:nyears)]))),
            sdCH4ha  = list(unname(unlist(resdf[r,65+c(1:nyears)]))),
            CH4_tyr = list(unname(unlist(resdf[r,107+c(1:nyears)]))),
            sdCH4  = list(unname(unlist(resdf[r,149+c(1:nyears)])))
          )
          restib<-rbind(restib,aoi_tib)
        }
        json_output<-jsonlite::toJSON(restib,
                                      auto_unbox = T,
                                      pretty = T)
        write(json_output, paste0("output/",intervention_type,".json"))
