# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /* 
# Filename: AnnualandPerennialCropland.R
# Author: Lisa Eash
# Date Started: 06/11/2024
# Last Edited: 06/24/2024
# Purpose: AFOLU GHG Calculations for Annual and Perennial Cropland Interventions
# 
# */

#### Required Inputs ####
    # Climate zone and soil type: These are contained in the geography_mask_area_info json
    #   file, which sources data from Dynamic World land use data,  data layers
    # User inputs: User inputs depend on subintervention category (tillage reduction,
    #   nutrient management, or livestock integration). 
    #   -For tillage reduction, only tillage type is required. 
    #   -For nutrient management, user must provide organic matter management (low, med,
    #     high, high with organic amendments) and synthetic N application rate and N%. If  
    #     the addition of organic amendments is part of the improved nutrient management,
    #     they may also enter organic amendment application rate and N%.
    #   -For livestock incorporation, user must provide organic matter management (low, med,
    #     high, high with organic amendments), livestock type, and stocking rate (head/ha).  

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
          # FPRP: Amount of N from grazing livestock waste (kg N/ha/yr)
          # EFDir: Emissions factor for direct N2O emissions from N inputs (kg N2O-N/kg N)
          # EF_PRP: Emissions factor for direct N2O emissions from grazing livestock waste (kg N2O-N/kg N)
          # FRAC_GASF: Fraction of synthetic fertilizer that volatilizes (kg NH3-N + NOx-N)/kg N
          # FRAC_GASM: Fraction of organic N that volatilizes (kg NH3-N + NOx-N)/kg N
          # EF_vol: Emissions factor for indirect N2O emissions from volatilization kg N2O-N/(kg NH3-N + NOx-N)
          # FRAC_LEACH: Fraction of N inputs that is lost to leaching (kg N loss/kg N)
          # EF_leach: Emissions factor for indirect N2O emissions from volatilization (kg N2O-N/kg N loss)
          # Nex: Nitrogen excretion rate from livestock (kg N/1000 kg animal mass/day)
      #CH4
          # EF_ch4: Emissions factor for CH4 from enteric fermentation from livestock incorporation (kg CH4/head/yr)

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
            j <- jsonlite::fromJSON('croplandex.json', flatten=TRUE)
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
            df$ef_direct_n2o_sd=(df$ef_direct_n2o_uncertainty_high-df$ef_direct_n2o_uncertainty_low)/3.92
            df$ef_leaching_sd=(df$ef_leaching_uncertainty_high-df$ef_leaching_uncertainty_low)/3.92
            df$frac_org_fert_sd=(df$frac_org_fert_uncertainty_high-df$frac_org_fert_uncertainty_low)/3.92
            df$frac_leach_sd=(df$frac_leach_uncertainty_high-df$frac_leach_uncertainty_low)/3.92
            df$frac_syn_fert_sd=(df$frac_syn_fert_uncertainty_high-df$frac_syn_fert_uncertainty_low)/3.92
            df$ef_prp_sd=(df$ef_prp_uncertainty_high-df$ef_prp_uncertainty_low)/3.92
            df$ef_vol_sd=(df$ef_vol_uncertainty_high-df$ef_vol_uncertainty_low)/3.92
            #Uncertainty presented as 95% CI as percentage of mean
            df$SOC_ref_tonnes_C_ha_sd=(df$high_activity_clay_soils_HAC_error_positive/100*df$SOC_ref_tonnes_C_ha)/1.96
            #Uncertainty presented as 2 sd as percentage of mean
            df$FLU_sd=(df$FLU_error_positive/100*df$FLU)/2
            df$FMG_sd=(df$FMG_error_positive/100*df$FMG)/2
            df$FI_sd=(df$FI_error_positive/100*df$FI)/2
            #Remove other uncertainty columns 
            df<-df[,-grep("uncertainty|error",names(df))]


####GHG Calculations####
    GHGcalc<-function(i,df,nx #number of MC iterations
    ){
      temp_bau<-df[df$aoi_id==i & df$scenario=="business-as-usual",]
      temp_int<-df[df$aoi_id==i & df$scenario=="intervention",]
      results.unc <- matrix(0, nx, 6)
      str(results.unc)
      results.unc[, 1] <- i
      results.unc[, 2] <-temp_bau$area
      results.unc[, 3] <- c(1:nx)
      colnames(results.unc) <- c('AOI','Area','Rep','SOC','N2O','CH4')
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
                if(grepl("livestock incorporation", interv_sub)){
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
                results.unc[m,5]<-(dirN+volN+leachN)/1000*44/28 #sum and convert to tN2O
        ####CH4 emissions####
            #emissions from enteric fermentation (if livestock incorporation)
            if(grepl("livestock incorporation", interv_sub)){
              EF_CH4<-rnorm(1, temp_int$ef_ch4, 
                            temp_int$ef_ch4*.50) #In Section 10.3.4 IPCC (2006), uncertainty range is 30-50% for livestock parameters
              CH4_live<-(temp_int$stocking_rate-temp_bau$stocking_rate)*EF_CH4/1000
            } else{ CH4_live<-0 }
        results.unc[m,6]<-CH4_live
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
        resdf <- matrix(0, length(AOIs), 14)
        resdf[, 1] <- AOIs
        colnames(resdf) <- c('AOI','Area','CO2_thayr','N2O_thayr','CH4_thayr',
                         'sdCO2ha','sdN2Oha','sdCH4ha',
                         'CO2_tyr','N2O_tyr','CH4_tyr',
                         'sdCO2','sdN2O','sdCH4')
        for(a in c(1:length(AOIs))){
          resdf[a, 2] <-mean(mc[[a]]$Area)
          resdf[a, 3] <- mean(mc[[a]]$SOC)
          resdf[a, 4] <- mean(mc[[a]]$N2O)
          resdf[a, 5] <- mean(mc[[a]]$CH4)
          resdf[a, 6] <- sd(mc[[a]]$SOC)
          resdf[a, 7] <- sd(mc[[a]]$N2O)
          resdf[a, 8] <- sd(mc[[a]]$CH4)
          resdf[a, 9] <- resdf[a,2]*resdf[a,3]
          resdf[a, 10] <- resdf[a,2]*resdf[a,4]
          resdf[a, 11] <- resdf[a,2]*resdf[a,5]
          resdf[a, 12] <- resdf[a,2]*resdf[a,6]
          resdf[a, 13] <- resdf[a,2]*resdf[a,7]
          resdf[a, 14] <- resdf[a,2]*resdf[a,8]
        }
        resdf<-as.data.frame(resdf[,-2])
    # convert to json output
        str(resdf)
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
