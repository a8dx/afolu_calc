# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: Livestock.R
# Author: Lisa Eash
# Date Started: 06/14/2024
# Last Edited: 06/29/2024
# Purpose: AFOLU GHG Calculations for Coastal Wetlands Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: The user must first select livestock included in their project from the following livestock categories:
    #   -Cattle:"cattle_dairy","cattle_male","cattle_female","cattle_calves"
    #   -Buffalo:"buffalo"
    #   -Wool: "goat_male","goat_female","goat_lamb","sheep_female","sheep_castrate","sheep_intact","sheep_lamb",
    #   -Swine: "swine_finishing","swine_breeding",
    #   -Poultry:"chicken","duck","turkey"
    #   For each livestock category (cattle, buffalo, wool, swine poultry), they must provide:
    #       -The percentage of manure managed in each animal waste management system (AWMS),
    #           -Pasture range and paddock
    #           -Solid storage
    #           -dry lot
    #           -uncovered anaerobic lagoon
    #           -anaerobic digester (type specified)
    #           -compost (type specified)
    #           -aerobic treatment
    #   For each livestock subcategory, they must provide:
    #       -Average live bodyweight (kg/head)
    #       -Number of head
    #       -Productivity type (high, low)
    #   Other user inputs depend on livestock category and include the following:
    #       -Wool:
    #           -feeding_situation: housed ewes, grazing flat pasture, grazing hilly pasture,
    #             housed fattening lambs, lowland goats, hill and mountain goats
    #           -weaning_weight: average weight at weaning (kg)
    #           -weight_slaughter: average weight at slaughter or 1 year (kg)
    #           -wool: are sheep/goats wool-producing (yes/no)
    #       -Cattle:
    #           -milk_production: average fat corrected milk production (kg/day) (dairy cattle only)
    #           -feed_type: High grain diet >90% (steam-flaked corn or other), high, moderate or low quality forage


#### Parameters and Paths ####
  # Tier 2 calculations are used for cattle and wool. Tier 1 is used for buffalo, swine, and poultry.
  # Calculations requires the following parameters. Where indicated, associated uncertainty is also considered:
      #Energy Intake- only for cattle and wool; swine, buffalo and poultry are calculated using Tier 1 emission factors
          # NEmf: estimated dietary net energy concentration of feed (MJ/kg DM; only for cattle-calves)
          # Wool:
              # Cfi: coefficient for net energy for maintenance (MJ/day/kg)
              # Ca: coefficient for net energy for activity based on animal's feeding situation (MJ/day/kg)
              # a_constant: constant a for net energy for growth
              # b_constant: constant b for net energy for growth
              # Cp: Coefficient for net energy for pregnancy
              # DE: Digestibility of feed, fraction of gross energy intake (%). Uncertainty 
                  # selected randomly from range of values in Table 10.2
      #CH4 from enteric fermentation
          # Ym: methane conversion factor (MJ/kg CH4; cattle and wool only). Uncertainty calculated from 
              #IPCC confidence intervals provided in table.
          # CH4_EF: emissions factor for CH4 from enteric fermentation (kg CH4/head/yr; swine and buffalo only). Uncertainty
              #calculated from a confidence interval of 40% as suggested in text of Ch. 10.
      #CH4 from manure management
          # VS_rate: Volatile solid excretion rate (kg VS/1000 kg animal mass/day; buffalo, swine and poultry only). Uncertainty
              #calculated from a confidence interval of 40% as suggested in text of Ch. 10.
          # DE: Digestibility of feed, fraction of gross energy intake (%; cattle only)
          # MCF: methane conversion factor for each AWMS, livestock, and climate combination (%). Uncertainty
              #calculated from a confidence interval of 30% as suggested in text of Ch. 10.
          # B0: max. methane producing capacity for manure produced by livestock category (m3 CH4/kg VS excreted). Uncertainty
              #calculated from a confidence interval of 40% as suggested in text of Ch. 10.
      #N2O from manure management
          # Nret: fraction of daily N intake retained by animal (cattle and wool only)
          # Crude protein: protein content of feed (%; cattle only)
          # Nex: annual N excretion rates (kg N/animal/yr; buffalo, swine and poultry only)
          # awms_EF_N2O: Direct N2O emissions factor for each AWMS (kg N2O-N/kg N). Uncertainty
              #calculated from a confidence interval of 100% as suggested in text of Ch. 10, except for the
              #case of pasture management which has reported 95% confidence interval in Table 11.1.
          # EF_N_vol: Indirect N2O emissions factor due to volatilization kg N2O-N/(kg NH3-N + NOx-N). Uncertainty
              #calculated from IPCC confidence intervals provided in table.
          # EF_N_leach: Indirect N2O emissions factor due to leaching (kg N2O-N/kg N loss). Uncertainty
              #calculated from IPCC confidence intervals provided in table.
          # frac_gas_awms: Fraction of N lost via volatilization for each AWMS (kg NH3-N + NOx-N)/kg N. Uncertainty
              # calculated from IPCC confidence intervals provided in tables.
          # frac_leach_awms: Fraction of N lost via leaching for each AWMS (kg N loss/kg N). 


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
    cattle<-c("cattle_dairy","cattle_male","cattle_female","cattle_calves")
    buffalo<-c("buffalo")
    wool<-c("goat_male","goat_female","goat_lamb","sheep_female",
            "sheep_castrate","sheep_intact","sheep_lamb")
    swine<-c("swine_finishing","swine_breeding")
    poultry<-c("chicken","duck","turkey")
    AWMSlist<-c("uncovered_anaerobic_lagoon",
                "solid_storage","dry_lot","pasture_range_paddock",
                "aerobic_treatment","compost",
                "anaerobic_digester")
    AWMSlistmod<-AWMSlist[!AWMSlist=="pasture_range_paddock"]
    AWMSlistmod2<-AWMSlist[!AWMSlist=="anaerobic_digester"]
    # Extract data from json input
        #Intervention parameters
            j <- jsonlite::fromJSON('livestock_ex.json', flatten=TRUE)
            interv_sub<- j$intervention_sub_category 
            livestock_subcat<- j$livestock_subcat 
            livestock_subcat<-unlist(strsplit(livestock_subcat,","))
            livestock_cat<- j$livestock_cat
            livestock_cat<-unlist(strsplit(livestock_cat,","))
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
        #Calculate mean for frac_gas for anaerobic digester (not provided in table)
            for(l_cat in livestock_cat){
              df$x<-(df[,names(df)==paste0(l_cat,"_frac_gas_anaerobic_digester_uncertainty_high")]+
                       df[,names(df)==paste0(l_cat,"_frac_gas_anaerobic_digester_uncertainty_low")])/2
              names(df)[names(df)=="x"]<-paste0(l_cat,"_frac_gas_anaerobic_digester")
            }

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
            #Remove other uncertainty columns 
            df<-df[,-grep("uncertainty",names(df))]
        #find mean for range values
            range<-grep("low",colnames(df))
            range_names<-colnames(df)[range]
            range_cols<-gsub("_low*.","\\1",range_names)
            for(b in range_cols){
              df$x<-(df[,names(df)==paste0(b,"_high")]+
                       df[,names(df)==paste0(b,"_low")])/2
              names(df)[names(df)=="x"]<-paste0(b)
            }
####GHG Calculations####
    GHGcalc<-function(i,df,nx #number of MC iterations
                      ){
      #Create results matrix
          temp<-df[df$aoi_id==i & df$scenario=="business-as-usual",]
          results.unc <- matrix(0, nx, 12)
          results.unc[, 1] <- i
          results.unc[, 2] <-temp$area
          results.unc[, 3] <- c(1:nx)
          colnames(results.unc) <- c('AOI','Area','Rep',
                                     'CO2_bau',
                                     'N2O_bau','CH4_bau',
                                     'CO2_int',
                                     'N2O_int','CH4_int',
                                     'CO2',
                                     'N2O','CH4')

      for(m in c(1:nx)){
        #Select random parameters for each monte carlo iteration
            #Common parameters
                FRAC_GASM<- rnorm(1, temp$FRAC_GASM, temp$FRAC_GASM_sd)
                FRAC_LEACH<- rnorm(1, temp$FRAC_LEACH, temp$FRAC_LEACH_sd)
                EF_LEACH<- rnorm(1, temp$EF_N_leach, temp$EF_N_leach_sd)
                EF_VOL<- rnorm(1, temp$EF_N_vol, temp$EF_N_vol_sd)
                Efdir_N2O_bs <-rnorm(1, temp$EFdir_N2O_bovine_swine, temp$EFdir_N2O_bovine_swine_sd)
                Efdir_N2O_wool <-rnorm(1, temp$EFdir_N2O_wool, temp$EFdir_N2O_wool_sd)
            #Parameters for each awms
                awms_par<-matrix(0,1,length(AWMSlistmod))
                colnames(awms_par)<-AWMSlistmod
                for(p in c(1:length(AWMSlistmod))){
                  a<-AWMSlistmod[p]
                  awms_par[,p]<-rnorm(1, temp[,names(temp)==paste0(a,"_EF_N2O")],
                                      temp[,names(temp)==paste0(a,"_EF_N2O")]/1.96)
                }
                awms_par<-as.data.frame(awms_par)
                awms_par <- type.convert(awms_par, as.is = T)
            #Parameters for each livestock subcat
                ls_par<-matrix(0,7,18)
                ls_par[,1]<-livestock_subcat
                colnames(ls_par)<-c("livestock_subcat","DE",paste0(AWMSlistmod,"_frac_gas"),
                                    "Ym","CH4_EF","VS_rate",paste0(AWMSlistmod2,"_MCF"),
                                    "b0")
                for(n in c(1:length(livestock_subcat))){
                  l<-livestock_subcat[n]
                  l_cat<-ifelse(l %in% wool, "wool",
                                ifelse(l %in% poultry,"poultry",
                                       ifelse(l %in% swine, "swine",
                                              ifelse(l=="cattle_dairy","cattle_dairy",
                                                     ifelse(l=="buffalo","buffalo",
                                                     "cattle_nondairy")))))
                  ls_par[n,2]<-ifelse(l_cat %in% c("wool","cattle_dairy","cattle_nondairy"),
                                      runif(1,temp[,names(temp)==paste0(l,"_DE_low")],
                                            temp[,names(temp)==paste0(l,"_DE_high")]),0)
                  for(p in c(1:length(AWMSlistmod))){
                    a<-AWMSlistmod[p]
                    ls_par[n,2+p]<-rnorm(1,temp[,names(temp)==paste0(l_cat,"_frac_gas_",a)],
                                         temp[,names(temp)==paste0(l_cat,"_frac_gas_",a,"_sd")])
                  }
                  ls_par[n,9]<-ifelse(l_cat %in% c("cattle_dairy","cattle_nondairy"),
                                      rnorm(1, temp[,names(temp)==paste0(l,"_Ym")],
                                            temp[,names(temp)==paste0(l,"_Ym")]*.20/1.96),0)
                  ls_par[n,9]<-ifelse(l_cat %in% c("wool"),
                                      rnorm(1, temp[,names(temp)==paste0(l,"_Ym")],
                                            temp[,names(temp)==paste0(l,"_Ym_sd")]/1.96),0)
                  ls_par[n,10]<-ifelse(l_cat %in% c("swine","buffalo"),
                                      rnorm(1, temp[,names(temp)==paste0(l_cat,"_CH4_ef_ent")],
                                            temp[,names(temp)==paste0(l_cat,"_CH4_ef_ent")]*.40/1.96),0)
                  ls_par[n,11]<-ifelse(l_cat %in% c("swine","poultry","buffalo"),
                                       rnorm(1, temp[,names(temp)==paste0(l,"_VS_rate")],
                                             temp[,names(temp)==paste0(l,"_VS_rate")]*.40/1.96),0)
                  for(p in c(1:length(AWMSlistmod2))){
                    a<-AWMSlistmod2[p]
                    ls_par[n,11+p]<-rnorm(1,temp[,names(temp)==paste0(a,"_CH4_MCF")],
                                          temp[,names(temp)==paste0(a,"_CH4_MCF")]*.30/1.96)
                  }
                  ls_par[n,18]<-rnorm(1, temp[,names(temp)==paste0(l_cat,"_b0_CH4")],
                                      temp[,names(temp)==paste0(l_cat,"_b0_CH4")]*.40/1.96)
                }
                ls_par<-as.data.frame(ls_par)
                ls_par <- type.convert(ls_par, as.is = T)
      #Calculations for each scenario
          for(z in c(1:length(scen))){
          s<-scen[z]
          temp<-df[df$aoi_id==i & df$scenario==s,]
          ls_res<-matrix(0,length(livestock_subcat),3)
          ls_res[,1]<-livestock_subcat
          colnames(ls_res)<-c("livestock_subcat","N2O","CH4")
          for(y in c(1:length(livestock_subcat))){
            l<-livestock_subcat[y]
            l_cat<-ifelse(l %in% wool, "wool",
                          ifelse(l %in% poultry,"poultry",
                                 ifelse(l %in% swine, "swine",
                                        ifelse(l=="cattle_dairy","cattle_dairy",
                                               ifelse(l=="buffalo","buffalo",
                                                      "cattle_nondairy")))))
          ####Energy Intake####
              #Cattle
                  if(l=="cattle_calves"){
                    DMI=(temp$cattle_calves_weight_kg^.75)*
                      ((0.0582*temp$cattle_calves_NEmf-
                          0.00266*(temp$cattle_calves_NEmf^2)-0.1128)/0.239*
                         temp$cattle_calves_NEmf)
                  }
                  if(l=="cattle_male"){
                    DMI=3.83+0.0143*temp$cattle_male_weight_kg*0.96
                  }
                  if(l=="cattle_female"){
                    DMI=3.184+0.01536*temp$cattle_male_weight_kg*0.96
                  }
                  if(l=="cattle_dairy"){
                    DMI=0.0185*temp$cattle_dairy_weight_kg+
                      0.305*temp$cattle_dairy_milk_production
                  }
                  if(l_cat %in% c("cattle_dairy","cattle_nondairy")){
                    GE=18.45*DMI
                  }
              #Sheep and goats
              if(l %in% wool){
                NEm=temp[,names(temp)==paste0(l,"_CFi")]*
                  (temp[,names(temp)==paste0(l,"_weight_kg")]^0.75) #Net energy for maintenance
                NEa=temp[,names(temp)==paste0(l,"_Ca")]*
                  temp[,names(temp)==paste0(l,"_weight_kg")] #Net energy for activity
                NEg=(temp[,names(temp)==paste0(l,"_weight_slaughter")]-
                       temp[,names(temp)==paste0(l,"_weaning_weight")])* #Net energy for growth
                        (temp[,names(temp)==paste0(l,"_a_constant")]+0.5*
                           temp[,names(temp)==paste0(l,"_b_constant")]*
                           (temp[,names(temp)==paste0(l,"_weight_slaughter")]+
                              temp[,names(temp)==paste0(l,"_weaning_weight")]))/365 
                if(l %in% c("sheep_female","goat_female")){
                  EVmilk<-ifelse(l=="sheep_female",4.6,3)
                  NEl=(5*(temp[,names(temp)==paste0(l,"_weight_slaughter")]-
                            temp[,names(temp)==paste0(l,"_weaning_weight")]))/365*EVmilk #net energy for lactation
                } else{ NEl=0 }
                if(temp[,names(temp)==paste0(l,"_wool")]=="yes"){
                  NEw=0.25 #net energy for wool production, default value sourced from Eq 10.12
                } else{NEw=0}
                if(l %in% c("sheep_female","goat_female")){
                  NEp=temp$Cp*NEm  #net energy for pregnancy
                } else{NEp=0}
                DE<-ls_par[ls_par$livestock_subcat==l,"DE"]
                REM=1.123-(.004092*DE)+(0.00001126*(DE^2))-25.4/DE
                REG=1.164-(.005*DE)+(0.00001308*(DE^2))-37.4/DE
                GE=((NEm+NEa+NEl+NEp)/REM+(NEg+NEw)/REG)/DE
                DMI=GE/18.45
              }
        ####CH4 from enteric fermentation#### 
            if(l %in% c(cattle,wool)){ #Tier 2 for cattle/wool
              EF_ent=GE*ls_par[ls_par$livestock_subcat==l,"Ym"]/100*365/55.65 #kg CH4/head/yr
            } else { if(l %in% c(swine, buffalo)){ #Tier 1 for swine and buffalo
              EF_ent= ls_par[ls_par$livestock_subcat==l,"CH4_EF"] 
              } else {EF_ent=0 #no enteric fermentation for poultry #kg CH4/head/yr
              } 
            }
            CH4_ef<-EF_ent*temp[,names(temp)==paste0(l,"_head")] #kg CH4/yr
        ####CH4 from manure management####
            UE=ifelse((l_cat %in% c("cattle_dairy","cattle_nondairy","buffalo") &&
                           temp[,names(temp)==paste0(l,"_feed_type")]=="high grain diet >90%" ) |
                        (l_cat=="wool" &&
                           temp[,names(temp)==paste0(l,"_feeding_situation")] %in% 
                           c("Housed ewes","Housed fattening lambs")),
                      0.02,0.04) #UE values pulled from Eq. 10.24
            ASH=0.06 #based on default value in Eq. 10.24
            if(l_cat %in% c("cattle_dairy","cattle_nondairy","wool")){ #Tier 2 methods for cattle and wool
              VS=(GE*(1-ls_par[ls_par$livestock_subcat==l,"DE"]/100)+ #kg VS/head/day
                    GE*UE)*((1-ASH)/18.45) 

            } else{ #Tier 1 for swine, buffalo and poultry
              VS=temp[,names(temp)==paste0(l,"_weight_kg")]/1000* #kg VS/head/day
                ls_par[ls_par$livestock_subcat==l,"VS_rate"] 
            }
            sumMCF<-vector()
            for(k in c(1:length(AWMSlistmod2))){ #For each manure management system
              awms<-AWMSlistmod2[k]
              MCF<-temp[,names(temp)==paste0(l_cat,"_mm_",awms,"_pct")]/100*
                ls_par[ls_par$livestock_subcat==l,names(ls_par)==paste0(awms,"_MCF")]
              sumMCF[k]<-MCF
            }
            MCFall<-sum(sumMCF)
            EF_mm=(VS*365)*ls_par[ls_par$livestock_subcat==l,"b0"]*.67*MCFall #kg CH4/head/yr
            CH4_mm=EF_mm*temp[,names(temp)==paste0(l,"_head")] #kg CH4/yr
        ####N2O from manure management####
            #Select crude protein %, calculate N excretion rate
                CP<-ifelse(l_cat %in% c("cattle_dairy","cattle_nondairy"),
                           temp[,names(temp)==paste0(l,"_crude_protein")],
                           15) #for wool, crude protein is sourced from study detailed in Annex 10B.3 (IPCC, 2019)
                if(l_cat %in% c("cattle_dairy","cattle_nondairy","wool")){
                  Nint=DMI*CP
                  Nex=Nint*(1-temp[,names(temp)==paste0(l,"_Nret")])*365 #kg N/head/yr
                  } else{ 
                  Nex=temp[,names(temp)==paste0(l,"_Nex")] #kg N/head/yr
                  }
            #Direct N2O emissions
                awms_n2o_sum<-vector()
                EF_N2O_pasture_range<-ifelse(l_cat %in% c("wool"),
                                             Efdir_N2O_wool,Efdir_N2O_bs)
                for(k in c(1:length(AWMSlistmod))){ #For each manure management system
                  awms<-AWMSlistmod[k]
                  awms_n2o<-temp[,names(temp)==paste0(l_cat,"_mm_",awms,"_pct")]/100*
                    awms_par[,names(awms_par)==awms]
                  awms_n2o_sum[k]<-awms_n2o
                }
                awms_n2o_sum[length(AWMSlistmod)+1]<-EF_N2O_pasture_range*
                  temp[,names(temp)==paste0(l_cat,"_mm_pasture_range_paddock_pct")]/100
                N2O_mm_dir<-sum(awms_n2o_sum)*Nex* 
                  temp[,names(temp)==paste0(l,"_head")]*44/28 #kg N2O/yr
            #Indirect N2O from Volatilization
                awms_nvol_sum<-vector()
                for(k in c(1:length(AWMSlistmod))){ #For each manure management system
                  awms<-AWMSlistmod[k]
                  awms_nvol<-temp[,names(temp)==paste0(l_cat,"_mm_",awms,"_pct")]/100*
                    ls_par[ls_par$livestock_subcat==l,names(ls_par)==paste0(awms,"_frac_gas")]
                  awms_nvol_sum[k]<-awms_nvol
                }
                awms_nvol_sum[length(AWMSlistmod)+1]<-FRAC_GASM*
                  temp[,names(temp)==paste0(l_cat,"_mm_pasture_range_paddock_pct")]/100
                N2O_mm_vol<-sum(awms_nvol_sum)*Nex*
                  temp[,names(temp)==paste0(l,"_head")]*EF_VOL*44/28 #kg N2O/yr
            #Indirect N2O from Leaching
                awms_nleach_sum<-vector()
                for(k in c(1:length(AWMSlistmod))){ #For each manure management system
                  awms<-AWMSlistmod[k]
                  awms_nleach<-temp[,names(temp)==paste0(l_cat,"_mm_",awms,"_pct")]/100*
                    temp[,names(temp)==paste0(l_cat,"_frac_leach_",awms)]
                  awms_nleach_sum[k]<-awms_nleach
                }
                awms_nleach_sum[length(AWMSlistmod)+1]<-FRAC_LEACH*
                  temp[,names(temp)==paste0(l_cat,"_mm_pasture_range_paddock_pct")]/100
                N2O_mm_leach<-sum(awms_nleach_sum)*Nex*
                  temp[,names(temp)==paste0(l,"_head")]*EF_LEACH*44/28 #kg N2O/yr
        ####Record results for each livestock type
          ls_res[y,2]<-(N2O_mm_dir+N2O_mm_vol+N2O_mm_leach)/1000 #t N2O/yr
          ls_res[y,3]<-(CH4_mm+CH4_ef) # t CH4/yr
          ls_res<-as.data.frame(ls_res)
          ls_res <- type.convert(ls_res, as.is = T)
          }
          results.unc[m,5+(z-1)*3]<-sum(ls_res$N2O)
          results.unc[m,6+(z-1)*3]<-sum(ls_res$CH4)
      }
      }
        results.unc[,10]<-results.unc[,7]-results.unc[,4]
        results.unc[,11]<-results.unc[,8]-results.unc[,5]
        results.unc[,12]<-results.unc[,9]-results.unc[,6]
        results.unc<-results.unc[,c(1:3,10:12)]
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
        resdf <- matrix(0, length(AOIs), 14)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area','CO2_thayr','N2O_thayr','CH4_thayr',
                             'sdCO2ha','sdN2Oha','sdCH4ha',
                             'CO2_tyr','N2O_tyr','CH4_tyr',
                             'sdCO2','sdN2O','sdCH4')
        for(a in c(1:length(AOIs))){
          resdf[a, 2] <-mean(mc[[a]]$Area)
          resdf[a, 9] <- mean(mc[[a]]$CO2)
          resdf[a, 10] <- mean(mc[[a]]$N2O)
          resdf[a, 11] <- mean(mc[[a]]$CH4)
          resdf[a, 12] <- sd(mc[[a]]$CO2)
          resdf[a, 13] <- sd(mc[[a]]$N2O)
          resdf[a, 14] <- sd(mc[[a]]$CH4)
          resdf[a, 3] <- resdf[a,9]/resdf[a,2]
          resdf[a, 4] <- resdf[a,10]/resdf[a,2]
          resdf[a, 5] <- resdf[a,11]/resdf[a,2]
          resdf[a, 6] <- resdf[a,12]/resdf[a,2]
          resdf[a, 7] <- resdf[a,13]/resdf[a,2]
          resdf[a, 8] <- resdf[a,14]/resdf[a,2]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        # convert to json output
        restib<-data.frame()
        for(r in c(1:nrow(resdf))){
          aoi_tib <-tibble(
            AOI = paste0("Sub_AOI-",r),
            CO2_thayr = list(rep(resdf[r,2],20)),
            sdCO2ha = list(rep(resdf[r,5],20)),
            CO2_tyr = list(rep(resdf[r,8],20)),
            sdCO2 = list(rep(resdf[r,11],20)),
            N2O_thayr = list(rep(resdf[r,3],20)),
            sdN2Oha = list(rep(resdf[r,6],20)),
            N2O_tyr = list(rep(resdf[r,9],20)),
            sdN2O = list(rep(resdf[r,12],20)),
            CH4_thayr = list(rep(resdf[r,4],20)),
            sdCH4ha  = list(rep(resdf[r,7],20)),
            CH4_tyr = list(rep(resdf[r,10],20)),
            sdCH4  = list(rep(resdf[r,13],20))
          )
          restib<-rbind(restib,aoi_tib)
        }
        json_output<-jsonlite::toJSON(restib,
                                      auto_unbox = T,
                                      pretty = T)
        write(json_output, paste0("output/",interv_sub,".json"))

