# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: LUCGrasslandtoCropland.R
# Author: Lisa Eash
# Date Started: 06/30/2024
# Last Edited: 06/30/2024
# Purpose: AFOLU GHG Calculations for Grassland Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: User inputs are required for both the grassland and cropland scenarios, and the user must 
    #   select whether the project is converting from grassland to cropland or cropland to grassland. 
    #   -Inputs for the grassland scenario are: 
    #       -grassland condition (nominally degraded or native, high intensity grazing, severely degraded, improved grassland)
    #       -number of improvements over initial scenario (0, 1 or 2)
    #       -synthetic N application rate (kg/ha) and N%
    #       -organic amendment application rate (kg/ha) and N%.
    #       -livestock type
    #       -stocking rate (head/ha)
    #       -average live weight of liveestock type (kg/head)
    #       -Frequency of burn (number of years)
    #   -Inputs for the cropland scenario are: 
    #       -Whether land conversion involved burning of biomass (yes, no)
    #       -crop type (annual or perennial)
    #       -organic matter management (low, medium, high, high with manure)
    #       -tillage type(full, reduced or no-till)
    #       -synthetic N application rate (kg/ha) and N%
    #       -organic amendment application rate (kg/ha) and N%.
    #       -livestock type
    #       -stocking rate (head/ha)
    #       -average live weight of liveestock type (kg/head)
#### Parameters and Paths ####
  # Calculation requires the following parameters for both scenarios and their associated uncertainty:
      #Biomass and dead OM C
          # crop_biomass: Amount of crop biomass lost or gained 1 year after conversion
          # grass_biomass: Amount of grassland biomass lost or gained 1 year after conversion
      #SOC
          # SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
          # FLU_grass: Land use emissions factor (always 1 for grasslands)
          # FMG: Management factor for both business as usual and intervention scenario 
          # FI: C input factor for both business as usual and intervention scenario 
      #N2O
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
    scen=c("cropland","")
    Cf_herbaceous <- 0.47 #IPCC (2006)

    # Extract data from json input
        #Intervention parameters
            j <- jsonlite::fromJSON('luc_grasslandcropland_ex.json', flatten=TRUE)
            interv_sub<- j$intervention_subcategory #intervention type
            common<-as.data.frame(j$scenarios$common)
            gdf<-as.data.frame(j$scenarios$grassland)
            gdf<-gdf[,!names(gdf)=="aoi_subregions"]
            gdf_aoi<-as.data.frame(j$scenarios$grassland$aoi_subregions)
            gdf_aoi$aoi_id<-c(1:nrow(gdf_aoi))
            gdf<-rbind(gdf,gdf)
            gdf<-cbind(gdf,gdf_aoi)
            gdf$scenario<-"grassland"
            cdf<-as.data.frame(j$scenarios$cropland)
            cdf<-cdf[,!names(cdf)=="aoi_subregions"]
            cdf_aoi<-as.data.frame(j$scenarios$cropland$aoi_subregions)
            cdf_aoi$aoi_id<-c(1:nrow(cdf_aoi))
            cdf<-rbind(cdf,cdf)
            cdf<-cbind(cdf,cdf_aoi)
            cdf$scenario<-"cropland"
            common<-rbind(common,common)
            cdf<-cbind(common,cdf)
            gdf<-cbind(common,gdf)
            cdf <- type.convert(cdf, as.is = T)
            gdf <- type.convert(gdf, as.is = T)
            AOIs<-unique(cdf$aoi_id)
        #convert uncertainty values to sd
          for(s in scen){
            if(s=="cropland"){
              df=cdf
            } else {df=gdf}
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
            if(s=="cropland"){
              cdf=df
            } else {gdf=df}
          }


####GHG Calculations####
    GHGcalc<-function(i,df,nx #number of MC iterations
    ){
      temp_g<-gdf[gdf$aoi_id==i,]
      temp_c<-cdf[cdf$aoi_id==i,]
      results.unc <- matrix(0, nx, 45)
      str(results.unc)
      results.unc[, 1] <- i
      results.unc[, 2] <-temp_g$area
      results.unc[, 3] <- c(1:nx)
      colnames(results.unc) <- c('AOI','Area','Rep','CO2_y1','CO2_y2',paste0('N2O_y',c(1:20)),
                                 paste0('CH4_y',c(1:20)))
      for(m in c(1:nx)){
        ####Biomass and Dead OM C#### - only year 1
            #Calculate loss of initial biomass in grassland
                grass_bio<-rnorm(1,temp_g$grass_biomass,temp_g$grass_biomass_sd)
                BioG<-Cf_herbaceous*grass_bio
            #Calculate gain in annual crop biomass 
                BioCrop<-rnorm(1,temp_c$crop_biomass,temp_c$crop_biomass_sd)
            #Calculate year 1 change in biomass
                if(interv_sub=="grassland to cropland"){
                  BioC_yr1=(BioG-BioCrop)*44/12
                }else{
                  BioC_yr1=(BioCrop-BioG)*44/12
                }
        ####SOC Stock Change####
                SOCREF<-rnorm(1, temp_g$SOC_ref_tonnes_C_ha, temp_g$SOC_ref_tonnes_C_ha_sd) #same for grass and crop
            #Calculate SOC of grassland
                FLUg<-temp_g$FLU_grass
                FMGg<-rnorm(1, temp_g$FMG_grass, temp_g$FMG_grass_sd)
                FIg<-rnorm(1, temp_g$FI_grass, temp_g$FI_grass_sd)
                SOCg<-SOCREF*FLUg*FMGg*FIg
            #Calculate SOC of cropland
                FLUc<-rnorm(1, temp_c$FLU_crop, temp_c$FLU_crop_sd)
                FMGc<-rnorm(1, temp_c$FMG_crop, temp_c$FMG_crop_sd)
                FIc<-rnorm(1, temp_c$FI_crop, temp_c$FI_crop_sd)
                SOCc<-SOCREF*FLUc*FMGc*FIc
            #Calculate difference
                if(interv_sub=="grassland to cropland"){
                  dSOC<-(SOCg-SOCc)/20
                }else{
                  dSOC<-(SOCc-SOCg)/20
                }
                results.unc[m,c(4)]<-dSOC*44/12+BioC_yr1
                results.unc[m,c(5)]<-dSOC*44/12
        ####N2O emissions####
            #Calculate change in N sources- N inputs in grassland vs. N inputs in cropland, assumed to be constant
                FSN_c<-temp_c$n_fertilizer_amount*temp_c$n_fertilizer_percent/100
                FSN_g<-temp_g$n_fertilizer_amount*temp_g$n_fertilizer_percent/100
                FON_c<-temp_c$org_amend_rate*temp_c$org_amend_npercent/100
                FON_g<-temp_g$org_amend_rate*temp_g$org_amend_npercent/100
                FSOM<-dSOC/10*1000 #FSOM 
                FPRP_c<-temp_c$live_weight/1000*temp_c$Nex*temp_c$stocking_rate*365
                FPRP_g<-temp_g$live_weight/1000*temp_g$Nex*temp_g$stocking_rate*365
                if(interv_sub=="grassland to cropland"){
                  FSN=FSN_c-FSN_g
                  FON=FSN_c-FSN_g
                  FPRP_g=-FPRP_g
                }else{
                  FSN=FSN_g-FSN_c
                  FON=FSN_g-FSN_c
                  FPRP_c=-FPRP_c
                }
            #Calculate N emissions
                EFdir<-rnorm(1,temp_g$ef_direct_n2o,temp_g$ef_direct_n2o_sd)
                FRAC_GASF<-rnorm(1,temp_g$frac_syn_fert,temp_g$frac_syn_fert_sd)
                FRAC_GASM<-rnorm(1,temp_g$frac_org_fert,temp_g$frac_org_fert_sd)
                EF_vol<-rnorm(1,temp_g$ef_vol,temp_g$ef_vol_sd)
                FRAC_LEACH<-rnorm(1, temp_g$frac_leach, temp_g$frac_leach_sd)
                EF_leach<-rnorm(1,temp_g$ef_leaching, temp_g$ef_leaching_sd)
                EF_PRP_g<-rnorm(1, temp_g$ef_prp, temp_g$ef_prp_sd)
                EF_PRP_c<-rnorm(1, temp_c$ef_prp, temp_c$ef_prp_sd)
                dirN<-(FSN+FON+FSOM)*EFdir + FPRP_c*EF_PRP_c + FPRP_g*EF_PRP_g #Direct emissions
                volN<-(FSN*FRAC_GASF+ (FON+FPRP_c+FPRP_g)*FRAC_GASM)*EF_vol #Indirect volatilization
                leachN<-(FSN+FON+FSOM+FPRP_c+FPRP_g)*FRAC_LEACH*EF_leach    #indirect leaching
                N2O<-(dirN+volN+leachN)/1000*44/28 #sum and convert to tN2O
                results.unc[m,c(6:25)]<-N2O
                #Add change in N2O due to fire management
                if(temp_g$fire_used=="True" | temp_c$conv_burning=="yes"){
                    fire_n2o_ef<-rnorm(1,temp_g$burning_n2o_ef_mean,temp_g$burning_n2o_ef_sd)
                    if(temp_c$conv_burning=="yes" && interv_sub=="cropland to grassland"){
                        CF_c<-rnorm(1,temp_c$combustion_factor_mean,
                                    temp_c$combustion_factor_mean*.4/1.96) #uncertainty value for ag. residues is not included in table,  assumed 40% uncertainty
                        MB_c<-BioCrop/Cf_herbaceous
                        fireN2O_c<-MB_c*CF_c*fire_n2o_ef/1000
                    }
                    if(temp_g$fire_used=="True"| (temp_c$conv_burning=="yes" && 
                                                  interv_sub=="grassland to cropland")){
                      CF_g<-rnorm(1,temp_g$combustion_factor_mean,temp_g$combustion_factor_sd)
                      MB_g<-rnorm(1,temp_g$fuel_biomass_mean,temp_g$fuel_biomass_sd)
                      fireN2O_g<-MB_g*CF_g*fire_n2o_ef/1000
                    }
                    if(temp_g$fire_used=="True" && interv_sub=="grassland to cropland"){
                      fire_per<-temp_g$fire_management_years
                      fire_yrs<-vector()
                      fire_yrs<-append(fire_yrs,1)
                      for (y in 1:19) { 
                        if (y%%fire_per == 0) {
                          fire_yrs<-append(fire_yrs,y+1)
                        }
                      }
                      results.unc[m,5+fire_yrs]<-N2O-fireN2O_g
                    }
                    if(temp_g$fire_used=="True" && interv_sub=="cropland to grassland"){
                      fire_per<-temp_g$fire_management_years
                      fire_yrs<-vector()
                      for (y in 1:19) { 
                        if (y%%fire_per == 0) {
                          fire_yrs<-append(fire_yrs,y+1)
                        }
                      }
                      results.unc[m,5+fire_yrs]<-results.unc[m,5+fire_yrs]+fireN2O_g
                    }
                    if(temp_c$conv_burning=="yes"){
                      if(interv_sub=="grassland to cropland"){
                        results.unc[m,6]<-results.unc[m,6]+fireN2O_g
                      }
                      if(interv_sub=="cropland to grassland"){
                        results.unc[m,6]<-results.unc[m,6]+fireN2O_c
                      }
                    }
                }
        ####CH4 emissions####
            #emissions from enteric fermentation- in grassland and then in cropland if livestock incorporated
              EF_CH4_c<-rnorm(1, temp_c$ef_ch4, 
                            temp_c$ef_ch4*.50/2) #In Section 10.3.4 IPCC (2006), uncertainty range (2 sd) is 30-50% for livestock parameters
              EF_CH4_g<-rnorm(1, temp_g$ef_ch4, 
                            temp_g$ef_ch4*.50/2)
              CH4_live_c<-(temp_c$stocking_rate)*EF_CH4_c/1000
              CH4_live_g<-(temp_g$stocking_rate)*EF_CH4_g/1000
              if(interv_sub=="grassland to cropland"){
                CH4_live_g=-CH4_live_g
              }else{
                CH4_live_c=-CH4_live_c
              }
              CH4_live=CH4_live_c+CH4_live_g
              results.unc[m,c(26:45)]<-CH4_live
            #emissions from fire management
            if(temp_g$fire_used=="True" | temp_c$conv_burning=="yes"){
              fire_ch4_ef<-rnorm(1,temp_g$burning_ch4_ef_mean,temp_g$burning_ch4_ef_sd)
              if(temp_c$conv_burning=="yes" && interv_sub=="cropland to grassland"){
                CF_c<-rnorm(1,temp_c$combustion_factor_mean,temp_c$combustion_factor_sd)
                MB_c<-BioCrop/Cf_herbaceous
                fireCH4_c<-MB_c*CF_c*fire_ch4_ef/1000
              }
              if(temp_g$fire_used=="True"| (temp_c$conv_burning=="yes" && 
                                            interv_sub=="grassland to cropland")){
                CF_g<-rnorm(1,temp_g$combustion_factor_mean,temp_g$combustion_factor_sd)
                MB_g<-rnorm(1,temp_g$fuel_biomass_mean,temp_g$fuel_biomass_sd)
                fireCH4_g<-MB_g*CF_g*fire_ch4_ef/1000
              }
              if(temp_g$fire_used=="True" && interv_sub=="grassland to cropland"){
                results.unc[m,25+fire_yrs]<-CH4_live-fireCH4_g
              }
              if(temp_g$fire_used=="True" && interv_sub=="cropland to grassland"){
                results.unc[m,25+fire_yrs]<-results.unc[m,25+fire_yrs]+fireCH4_g
              }
              if(temp_c$conv_burning=="yes"){
                if(interv_sub=="grassland to cropland"){
                  results.unc[m,26]<-results.unc[m,26]+fireCH4_g
                }
                if(interv_sub=="cropland to grassland"){
                  results.unc[m,26]<-results.unc[m,26]+fireCH4_c
                }
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
        resdf <- matrix(0, length(AOIs), 170)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area','CO2_thayr_y1','CO2_thayr_y2',
                             paste0('N2O_thayr_y',c(1:20)),
                             paste0('CH4_thayr_y',c(1:20)),
                             'sdCO2ha_y1','sdCO2ha_y2',paste0('sdN2Oha_y',c(1:20)),
                             paste0('sdCH4ha_y',c(1:20)),
                             'CO2_tyr_y1','CO2_tyr_y2',
                             paste0('N2O_tyr_y',c(1:20)),
                             paste0('CH4_tyr_y',c(1:20)),
                             'sdCO2_y1','sdCO2_y2',
                             paste0('sdN2O_y',c(1:20)),
                             paste0('sdCH4_y',c(1:20)))
        for(a in c(1:length(AOIs))){
          temp_res<-as.data.frame(mc[[a]])
          resdf[a, 2] <-mean(temp_res$Area)
          for(z in c(4:45)){
            resdf[a, z-1] <- mean(temp_res[,z])
            resdf[a, z+41] <- sd(temp_res[,z])
          }
        }
        for(k in c(3:86)){
          resdf[, 84+k] <- resdf[,2]*resdf[,k]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        restib<-data.frame()
        for(r in c(1:nrow(resdf))){
          aoi_tib <-tibble(
            AOI = paste0("Sub_AOI-",r),
            CO2_thayr = list(c(resdf[r,2],rep(resdf[r,3],19))),
            sdCO2ha = list(c(resdf[r,44],rep(resdf[r,45],19))),
            CO2_tyr = list(c(resdf[r,86],rep(resdf[r,87],19))),
            sdCO2 = list(c(resdf[r,128],rep(resdf[r,129],19))),
            N2O_thayr = list(unname(unlist(resdf[r,c(4:23)]))),
            sdN2Oha = list(unname(unlist(resdf[r,c(46:65)]))),
            N2O_tyr = list(unname(unlist(resdf[r,c(88:107)]))),
            sdN2O = list(unname(unlist(resdf[r,c(130:149)]))),
            CH4_thayr = list(unname(unlist(resdf[r,c(24:43)]))),
            sdCH4ha  = list(unname(unlist(resdf[r,c(66:85)]))),
            CH4_tyr = list(unname(unlist(resdf[r,c(108:127)]))),
            sdCH4  = list(unname(unlist(resdf[r,c(150:169)])))
          )
          restib<-rbind(restib,aoi_tib)
        }
        json_output<-jsonlite::toJSON(restib,
                                      auto_unbox = T,
                                      pretty = T)
        write(json_output, paste0("output/",interv_sub,".json"))
