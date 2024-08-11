# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: Forest_Protection_Deforestation_test.py
# Author: Barbara Bomfim
# Date Started: 07/11/2024
# Last Edited: 07/19/2024 (Anthony D'Agostino)
# Purpose: AFOLU GHG Calculations for Forest Protection from Deforestation (FP) Interventions
# Updates:
#      7/19/2024 - switch from L4A to L4B as AGB source
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data_FP from Dynamic World land use data_FP,  data_FP layers
# User inputs: User inputs depend on sub-intervention category (planted or natural forest; 
#  tillage; nutrient management; fire management). 
#   - Initial:
#           -land area: area under business-as-usual land use
#           -forest_type_deforestation: forest type prior to deforestation
#.          -forest_management_type: Natural, Plantation
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#           -fire_used: yes, no
#.  -Baseline data_FP:
#           -deforested_land_area: area under business-as-usual land use (deforestation)
#           -deforestation_rate_pre: Deforestation rate (in %/yr) before intervention, based on published data_FP (e.g., global forest watch)
#           -forest_type_deforestation: forest type prior to deforestation
#.          -forest_management_type: Natural, Plantation
#           -fire_used: yes, no
#           -deforested_for_fuelwood: Percentage of harvested wood used for fuelwood
#           -agroforestry_type: Alleycropping, Fallow, Hedgerow, Multistrata, Parkland, Shaded perennial, Silvoarable, Silvopasture
#           -monoculture_type: Oilpalm, Rubber, Tea
#           -settlement_type: Settlement, Urban green, Turfgrass, Cultivated soil
#           -reservoir_type: Default, Oligotrophic, Mesotrophic, Eutrophic, Hypereutrophic
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#           -illegal_logging_rate: volume timber over bark extracted illegally each year (m3 ha-1 yr-1)
#           -D: wood density (t/m3)
#  -Intervention requires:
#           -forest_area: forest area that is being protected
#           -forest_type: selection from dropdown list
#           -deforestation_rate_post: Deforestation rate (in %/yr) after intervention, based on published data_FP (e.g., global forest watch)
#           -fire_used: yes, no
#.          -fire_incidence_rate: default annual average fire incidence rate (0.1) but can be calculated as the average of the years of available data_FP. 
#               Calculate fraction of forest (by land cover data_FP) that had burned each year to get the annual fire incidence rate.
#           -bio_ef: gCO2 per kg of dry matter burnt
#           -D: wood density (t/m3)

#### Parameters and Paths ####
# Calculations requires the following parameters and their associated uncertainty:

##AGB
#AGBREF:AGB for the climate zone and soil type (tC/ha)

##BGB
#BGBREF:BGB for the climate zone and soil type (tC/ha)

##SOC
# SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
# FLU: Land use emissions factor (1 for planted forest, 1 for natural forest)
# FMG: Management factor for both business as usual scenario 
# FI: C input factor for both business as usual and intervention scenario 

##N2O
# burning_n2o_ef: Emissions factor for N2O from burning (g N2O/kg dry matter burnt)
# combustion_factor: combustion factor for fire
# fuel_biomass: amount of biomass available for burning (t dry matter/ha)
# FON: Amount of organic N applied annually (kg N/ha/yr)
# EFDir: Emissions factor for direct N2O emissions from N inputs (kg N2O-N/kg N)
# FRAC_GASM: Fraction of organic N that volatilizes (kg NH3-N + NOx-N)/kg N
# EF_vol: Emissions factor for indirect N2O emissions from volatilization kg N2O-N/(kg NH3-N + NOx-N)
# FRAC_LEACH: Fraction of N inputs that is lost to leaching (kg N loss/kg N)
# EF_leach: Emissions factor for indirect N2O emissions from volatilization (kg N2O-N/kg N loss)

##CH4
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
#     "N2O_tyr": [11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, 11.6191, #     11.6191, 11.6191, 11.6191, 11.6191],
#     "sdN2O": [17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206, 17.0206],
#     "CH4_thayr": [1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491, 1.1491],
#     "sdCH4ha": [0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849, 0.5849],
#     "CH4_tyr": [5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457, 5.7457],
#     "sdCH4": [2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247, 2.9247]
#   }
# ]

## Annual GHG impact GEDI data_FP and further calculations ###

## Annual GHG impact GEDI data_FP and further calculations ###

## Import libraries
import json
import os

import pandas as pd
import numpy as np
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor
import ee
import re
import math

from helpersFP import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock, sanitize_filename
#from helpersFP import calculate_historical_deforestation_rate, calculate_area_deforested, calculate_carbon_stock_loss, calculate_annual_baseline_emissions
#from helpersFP import calculate_ELE, calculate_Ldf_FP, calculate_incidental_damage


### STEP 1 - Estimate mean annual AGB (in Mg/ha) using GEDI data_FPset ####
# GEE path: https://code.earthengine.google.com/?scriptPath=users%2Fbabomfimf%2FAGB-GEDI-L4A%3Atest-3_AGB_annual_mean
## Anthony's GEE: for Persistent: `021_extract_SOC_and_AGB.js` in the `AFOLU` shared GEE Scripts library (for others: https://code.earthengine.google.com/44d6aaa5db4764b5b5f3825baf900d04)

#if this doesn't work authenticate from the command line  by running `earthengine authenticate`
ee.Authenticate()
# Initialize Earth Engine
# Use the name of GEE project here
ee.Initialize(project='ee-babomfimf3')

# Define sample polygon vertices
polygon_vertices = [
    [-55.85747628925944, -15.391799301778946],
    [-55.858920665940275, -15.411804954738642],
    [-55.83219969735529, -15.407813261205916],
    [-55.85747628925944, -15.391799301778946]
    ]
# Create an Earth Engine Polygon
polygon = ee.Geometry.Polygon(polygon_vertices)

# Access the GEDI L4A Monthly data_FPset
gedi_l4b = ee.Image('LARSE/GEDI/GEDI04_B_002')
              
# Filter the collection to your polygon
filtered_gedi = gedi_l4b.filterBounds(polygon)
# Filter agbd
agbd_image = filtered_gedi.select('MU')

# Compute the mean annual AGBD value for the entire polygon
average_agbd = agbd_image.reduceRegion(
    reducer=ee.Reducer.mean(),
    geometry=polygon,
    scale=1000  # Scale should match the resolution of the data_FPset
)
average_agbd_dict = average_agbd.getInfo()
print(f"{average_agbd_dict} Mg/ha")

print("\nArea of the polygon:")
print(f"{polygon.area().getInfo()/10000} ha") # to get area in hectares
maskedArea_ha = polygon.area().getInfo()/10000


## JSON data_FP input #####

# JSON file path - Get the current working directory
current_directory = os.getcwd()
print("Current working directory:", current_directory)
# Define the relative file path
file_name_FP = 'FP_test.json'
file_path_FP = os.path.join(current_directory, file_name_FP)
json_file_path_FP = './FP_test.json'

# Load data_FP from the JSON file
with open(file_path, 'r') as file:
        data = json.load(file)
data  

####STEP 2 #########
### Calculate mean annual AGC stock in tCO2e/ha####

def mean_annual_agc_stock(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    
    average_agbd_value = average_agbd_dict.get('agbd_mean', 0)
    area_converted = polygon.area().getInfo()/10000
    area_converted_yr = float(area_converted)
    
    CF = scenario_data["CF"]
    average_agbd_cstock = convert_to_c(average_agbd_value, CF)
    average_agbd_tco2e = convert_to_co2e(average_agbd_cstock)
    
    if log_level == 'debug':
        print(f"{area_converted} ha")
        print(f"Average Aboveground Biomass Density: {average_agbd_value} Mg/ha")
        print(f"Average Aboveground Carbon Stock: {average_agbd_tco2e} tCO2e/ha")
    
    subregion_results = []
    for subregion in scenario_data["aoi_subregions"]:
        subregion_area = subregion["area"]
        subregion_fraction = subregion_area / area_converted_yr
        
        subregion_agbd_tco2e = average_agbd_tco2e * subregion_fraction
        
        subregion_results.append({
            "aoi_id": subregion["aoi_id"],
            "area": subregion_area,
            "carbon_stock": subregion_agbd_tco2e
        })
        
        if log_level == 'debug':
            print(f"Subregion {subregion['aoi_id']}: {subregion_agbd_tco2e} tCO2e")
    
    return subregion_results
  
####STEP 3 ############
### Calculate mean annual Total C stock (sum AGC and BGC) ####

def mean_annual_tot_c_stock(subregion_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
    
    total_results = []
    
    for subregion in subregion_results:
        average_agbd_tco2e = subregion["carbon_stock"]
        
        average_bgbd_tco2e = convert_to_bgb(average_agbd_tco2e, ratio)
        average_total_tco2e = average_agbd_tco2e + average_bgbd_tco2e
        
        total_results.append({
            "aoi_id": subregion["aoi_id"],
            "area": subregion["area"],
            "aboveground_carbon_stock": average_agbd_tco2e,
            "belowground_carbon_stock": average_bgbd_tco2e,
            "total_carbon_stock": average_total_tco2e
        })
        
        if log_level == 'debug':
            print(f"Subregion {subregion['aoi_id']}:")
            print(f"  Aboveground Carbon Stock: {average_agbd_tco2e:.2f} tCO2e/ha")
            print(f"  Belowground Carbon Stock: {average_bgbd_tco2e:.2f} tCO2e/ha")
            print(f"  Total Carbon Stock: {average_total_tco2e:.2f} tCO2e/ha")
    
    return total_results

####STEP 4 #############
### BASELINE: Annual Total Biomass CO2 calculations ###

## Equations:
# historical_deforestation_rate = calculate_historical_deforestation_rate(deforestation_rate_pre, deforestation_rate_post)
# area_deforested_value = calculate_area_deforested(historical_deforestation_rate, forest_area)
# carbon_stock_loss = calculate_carbon_stock_loss(area_deforested_value, average_total_tco2e)
# annual_baseline_biomass_co2_emissions = calculate_annual_baseline_emissions(historical_deforestation_rate, carbon_stock_loss)


# Define equation to calculate Historical Deforestation Rate 
  # this will be provided from Global Forest Watch data
def calculate_historical_deforestation_rate(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        deforestation_pre = scenario_data.get("deforestation_rate_pre")
        deforestation_post = scenario_data.get("deforestation_rate_post")
    
    historical_def_rate = deforestation_pre - deforestation_post
    
    total_results.append({
        "hist_def_rate": historical_def_rate
    })
    
    if log_level == 'debug':
        # print(f"Subregion {subregion['aoi_id']}:")
        print(f"  Historical Deforestation Rate: {historical_def_rate:.2f} %/yr")
    
    return total_results

# Define equation to calculate Baseline Deforested Area
def calculate_deforested_area_yeari(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])

    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area_aoi = subregion["area"]
        historical_def_rate_aoi = subregion["hist_def_rate"]
    
    deforested_area_yeari = (historical_def_rate_aoi / 100) * forest_area_aoi
    
    total_results.append({
        "deforested_area_year_i": deforested_area_yeari
    })
    
    if log_level == 'debug':
        # print(f"Subregion {total_results['aoi_id']}:")
        print(f"  Deforested Area year i: {deforested_area_yeari:.2f} ha/yr")
    
    return total_results

# Define equation to calculate Baseline Total Biomass Carbon Stock Loss
def calculate_carbon_stock_loss_yeari(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    #ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        deforested_area_yeari = subregion["deforested_area_year_i"]
        average_total_tco2e = subregion["total_carbon_stock"]
 
        carbon_stock_loss = deforested_area_yeari * average_total_tco2e
        
        total_results.append({
            "total_carbon_stock_loss": carbon_stock_loss
        })
        
        if log_level == 'debug':
            print(f"Subregion {total_results['aoi_id']}:")
            print(f"  Carbon Stock Loss year i: {carbon_stock_loss:.2f} tCO2/yr")
    
    return total_results

# Define equation to calculate Baseline Total Biomass CO2 Emissions ##
def calculate_baseline_biomass_co2_emissions(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    #ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        historical_def_rate = subregion["hist_def_rate"]
        carbon_stock_loss_yeari = subregion["total_carbon_stock_loss"]
 
        annual_baseline_biomass_co2_emissions = historical_def_rate * carbon_stock_loss_yeari
        
        total_results.append({
            "annual_baseline_biomass_co2": annual_baseline_biomass_co2_emissions
        })
        
        if log_level == 'debug':
            print(f"Subregion {total_results['aoi_id']}:")
            print(f"  Annual Baseline Emissions year i: {annual_baseline_biomass_co2_emissions:.2f} tCO2/yr")
    
    return total_results

####### Baseline emissions calculations ###########

## Calculate all equations for baseline emissions calculation
def calculate_all(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    historical_def_rate = calculate_historical_deforestation_rate(deforestation_pre, deforestation_post, scenario=scenario, log_level=log_level)
    result["historical_def_rate"] = historical_def_rate

    deforested_area = calculate_deforested_area_yeari(historical_def_rate, forest_area_aoi, scenario=scenario, log_level=log_level)
    result["deforested_area"] = deforested_area

    carbon_stock_loss = calculate_carbon_stock_loss_yeari(deforested_area_yeari, average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["carbon_stock_loss"] = carbon_stock_loss

    annual_baseline_biomass_co2_emissions = calculate_baseline_biomass_co2_emissions(historical_def_rate, carbon_stock_loss, scenario=scenario, log_level=log_level)
    result["annual_baseline_biomass_co2_emissions"] = annual_baseline_biomass_co2_emissions
    if log_level == 'debug':
        print(f"Baseline Annual Biomass CO2 emissions: {annual_baseline_biomass_co2_emissions} tCO2e/yr")

    return result
result["annual_baseline_biomass_co2_emissions"]



####### Intervention Calculations: Avoided Deforestation Total Biomass CO2 Calculations #######

#Define equations needed to calculate Avoided Deforestation Emissions from Trees:
#Eq 1 - Actual deforested area and Forest area end of year n
def calculate_actual_area_deforested_year_n(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    forest_area = scenario_data["forest_area"]
    deforestation_post = scenario_data["deforestation_rate_post"]
    
    total_results = []
    
    forest_area_end_of_year_n = forest_area

    for year in range(1, n + 1):
        actual_area_deforested_year_n = (forest_area_end_of_year_n * deforestation_post) / 100
        forest_area_end_of_year_n -= actual_area_deforested_year_n
    
    total_results.append({
            "actual_area_deforested_year_n": actual_area_deforested_year_n,
            "forest_area_end_of_year_n": forest_area_end_of_year_n
    })
        
    if log_level == 'debug':
        print(f"Subregion {total_results['aoi_id']}:")
        print(f"actual_area_deforested_year_n: {actual_area_deforested_year_n:.2f} ha/yrr")
        print(f"forest_area_end_of_year_n: {forest_area_end_of_year_n:.2f} ha/yrr")

    return total_results

#Eq 2 - area of avoided deforestation in year n
def calculate_area_avoided_deforestation_year_n(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    forest_area = scenario_data["forest_area"]
    deforestation_post = scenario_data["deforestation_rate_post"]
    
    total_results = []
    
    forest_area_each_year = [forest_area]
    
    for year in range(1, years + 1):
        actual_deforested_area = forest_area_each_year[-1] * deforestation_post / 100
        forest_area_end_of_year = forest_area_each_year[-1] - actual_deforested_area
        forest_area_each_year.append(forest_area_end_of_year)

    total_results.append({
            "actual_deforested_area": actual_deforested_area,
            "forest_area_end_of_year": forest_area_end_of_year
    })
        
    if log_level == 'debug':
        print(f"Subregion {total_results['aoi_id']}:")
        print(f"forest_area_end_of_year: {forest_area_end_of_year:.2f} ha/yrr")
    
    return total_results

#eq 4 - Avoided Deforestation Emissions from Trees
def calculate_avoided_emissions_trees_each_year(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    forest_area = scenario_data["forest_area"]
    deforestation_post = scenario_data["deforestation_rate_post"]

    total_results = []
    
    forest_area_each_year = [forest_area]

    for year in range(1, years + 1):
        actual_deforested_area = forest_area_each_year[-1] * deforestation_post / 100
        forest_area_end_of_year = forest_area_each_year[-1] - actual_deforested_area
        forest_area_each_year.append(forest_area_end_of_year)

        avoided_deforestation_area = forest_area_each_year[-2] * deforestation_post / 100
        avoided_emissions_trees = avoided_deforestation_area * average_total_tco2e * (44 / 12)
        avoided_emissions_each_year.append(avoided_emissions_trees)

        total_results.append({
            "avoided_emissions_trees": avoided_emissions_trees,
            "forest_area_end_of_year": forest_area_end_of_year
         })
        
        if log_level == 'debug':
            print(f"Subregion {total_results['aoi_id']}:")
            print(f"avoided_emissions_trees: {avoided_emissions_trees:.2f} tCO2/yr")
     
        return total_results

######Calculate all - equations for baseline and intervention emissions calculation
def calculate_all(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    historical_def_rate = calculate_historical_deforestation_rate(deforestation_pre, deforestation_post, scenario=scenario, log_level=log_level)
    result["historical_def_rate"] = historical_def_rate

    deforested_area = calculate_deforested_area_yeari(historical_def_rate, forest_area_aoi, scenario=scenario, log_level=log_level)
    result["deforested_area"] = deforested_area

    carbon_stock_loss = calculate_carbon_stock_loss_yeari(deforested_area_yeari, average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["carbon_stock_loss"] = carbon_stock_loss

    annual_baseline_biomass_co2_emissions = calculate_baseline_biomass_co2_emissions(historical_def_rate, carbon_stock_loss, scenario=scenario, log_level=log_level)
    result["annual_baseline_biomass_co2_emissions"] = annual_baseline_biomass_co2_emissions
    if log_level == 'debug':
        print(f"Baseline Annual Biomass CO2 emissions: {annual_baseline_biomass_co2_emissions} tCO2e/yr")

    actual_area_deforested_each_year = calculate_actual_area_deforested_year_n(scenario=scenario, log_level=log_level)
    result["actual_area_deforested_each_year"] = actual_area_deforested_each_year
    
    forest_area_each_year = calculate_forest_area_end_of_year_n(scenario=scenario, n=20)
    result["forest_area_each_year"] = forest_area_each_year
    
    area_avoided_deforestation_year_n = calculate_area_avoided_deforestation_year_n(scenario=scenario, log_level=log_level)
    result["area_avoided_deforestation_year_n"] = area_avoided_deforestation_year_n
    
    yearly_avoided_emissions_trees = calculate_avoided_emissions_trees_each_year(scenario=scenario, log_level=log_level)
    result["yearly_avoided_emissions_trees"] = yearly_avoided_emissions_trees
    
    return result
result["yearly_avoided_emissions_trees"]

# ######## FOREST PROTECTION INTERVENTION CO2 IMPACT #########

## Calculate difference between intervention and business_as_usual
biomass = {}
biomass["business_as_usual"] = calculate_all(scenario="business_as_usual")
biomass["intervention"] = calculate_all(scenario="intervention")

biomass_co2_result = {}

for inter_item, bau_item in zip(biomass["intervention"]["yearly_avoided_emissions_trees"], 
                                biomass["business_as_usual"]["yearly_avoided_emissions_trees"]):
    if inter_item['aoi_id'] != bau_item['aoi_id']:
        raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
    
    aoi_id = inter_item['aoi_id']
    difference = inter_item['annual_change'] - bau_item['annual_change']
    
    biomass_co2_result[aoi_id] = difference# inter_result = calculate_all(scenario="intervention")

print(biomass_co2_result)

biomass_co2_result_error_positive = 10  # this is a percentage

biomass_co2_sd = {}

for aoi_id, difference in biomass_co2_result.items():
    # Calculate SD for each aoi_id
    sd = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96
    biomass_co2_sd[aoi_id] = sd

######## SOIL GHG CALCULATIONS #########

## Getting JSON data ready to be calculated ##
def load_json_data(file_path):
    with open(file_path, 'r') as f:
        json_data = json.load(f)
    
    interv_sub = json_data['intervention_subcategory']
    common = pd.DataFrame(json_data['scenarios']['common'])
    
    # Business as usual scenario
    bau = pd.DataFrame(json_data['scenarios']['business_as_usual'])
    bau = bau.drop('aoi_subregions', axis=1)
    bau_aoi = pd.DataFrame(json_data['scenarios']['business_as_usual'][0]['aoi_subregions'])
    # NOTE: these are not required with the new aoi_id proposed
    # bau_aoi['aoi_id'] = range(1, len(bau_aoi) + 1)
    bau = pd.concat([bau] * len(bau_aoi))  # Repeat bau data for each aoi
    bau = pd.concat([bau.reset_index(drop=True), bau_aoi.reset_index(drop=True)], axis=1)
    bau['scenario'] = 'business-as-usual'
    
    # Intervention scenario
    int_data = pd.DataFrame(json_data['scenarios']['intervention'])
    int_data = int_data.drop('aoi_subregions', axis=1)
    int_aoi = pd.DataFrame(json_data['scenarios']['intervention'][0]['aoi_subregions'])
    # NOTE: these are not required with the new aoi_id proposed
    # int_aoi['aoi_id'] = range(1, len(int_aoi) + 1)
    int_data = pd.concat([int_data] * len(int_aoi))  # Repeat int_data for each aoi
    int_data = pd.concat([int_data.reset_index(drop=True), int_aoi.reset_index(drop=True)], axis=1)
    int_data['scenario'] = 'intervention'
    
    # Combine all data
    df = pd.concat([bau, int_data])
    common = pd.concat([common] * len(df))
    df = pd.concat([common.reset_index(drop=True), df.reset_index(drop=True)], axis=1)
    
    # Convert data types
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='ignore')
    
    AOIs = df['aoi_id'].unique()
    
    # Convert uncertainty values to sd
    unc_low_cols = [col for col in df.columns if 'uncertainty_low' in col]
    unc_cols = [col.replace('_uncertainty_low', '') for col in unc_low_cols]
    
    for col in unc_cols:
        df[f'{col}_sd'] = (df[f'{col}_uncertainty_high'] - df[f'{col}_uncertainty_low']) / 3.92

    print(df)
    
    df = df.drop(columns=[col for col in df.columns if 'uncertainty' in col])
    
    # Handle specific uncertainty calculations
    df['SOC_ref_tonnes_C_ha_sd'] = (df['high_activity_clay_soils_HAC_error_positive'] / 100 * df['SOC_ref_tonnes_C_ha']) / 1.96
    df['FMG_sd'] = (df['FMG_error_positive'] / 100 * df['FMG']) / 2
    df['FI_sd'] = (df['FI_error_positive'] / 100 * df['FI']) / 2
    
    # Remove other uncertainty columns
    df = df.drop(columns=[col for col in df.columns if 'error' in col])
    
    return json_data, df, AOIs, interv_sub

###############################################################

def GHGcalc(aoi_id, df, nx, intervention_subcategory, biomass_co2_result):
    
    temp_bau = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "business-as-usual")]
    temp_int = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "intervention")]
    # Create a structured array that can hold both strings and floats
    results_unc = {}
    
    results_unc['AOI'] = aoi_id
    results_unc['Area'] = temp_bau['area'].values[0]
    results_unc['Rep'] = np.arange(1, nx + 1)
    results_unc['SOC'] = []
    results_unc['totalC'] = []
    results_unc[f'N2O_{aoi_id}'] = []
    results_unc[f'CH4_{aoi_id}'] = []
     
    
    for m in range(nx):
        # SOC Stock Change
        SOCREF = np.random.normal(temp_bau['SOC_ref_tonnes_C_ha'].values[0], temp_bau['SOC_ref_tonnes_C_ha_sd'].values[0])
        FLUbau = temp_bau['FLU'].values[0]
        FMGbau = np.random.normal(temp_bau['FMG'].values[0], temp_bau['FMG_sd'].values[0])
        FIbau = np.random.normal(temp_bau['FI'].values[0], temp_bau['FI_sd'].values[0])
        FLUint = temp_int['FLU'].values[0]
        FMGint = np.random.normal(temp_int['FMG'].values[0], temp_int['FMG_sd'].values[0])
        FIint = np.random.normal(temp_int['FI'].values[0], temp_int['FI_sd'].values[0])
        SOCbau = SOCREF * FLUbau * FMGbau * FIbau
        SOCint = SOCREF * FLUint * FMGint * FIint
        dSOC = SOCbau - SOCint
        results_unc["SOC"].append(dSOC)

        # here we add the biomass specific to the aoi_id
        # NOTE: this number is too large check please
        results_unc["totalC"].append(dSOC * 44/12 + biomass_co2_result[aoi_id]) 

        # N2O emissions
        # Calculate change in N sources
        FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
               temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100) if "nutrient management" in intervention_subcategory else 0
        FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
               temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100) if "nutrient management" in intervention_subcategory else 0
        FSOM = dSOC / 10 * 1000  # FSOM
        if "change in livestock type or stocking rate" in intervention_subcategory:
            FPRP = (temp_int['live_weight'].values[0] * temp_int['Nex'].values[0] * temp_int['stocking_rate'].values[0] -
                    temp_bau['live_weight'].values[0] * temp_bau['Nex'].values[0] * temp_bau['stocking_rate'].values[0]) / 1000 * 365
        else:
            FPRP = 0

        # Calculate N emissions
        EFdir = np.random.normal(temp_bau['ef_direct_n2o'].values[0], temp_bau['ef_direct_n2o_sd'].values[0])
        FRAC_GASF = np.random.normal(temp_bau['frac_syn_fert'].values[0], temp_bau['frac_syn_fert_sd'].values[0])
        FRAC_GASM = np.random.normal(temp_bau['frac_org_fert'].values[0], temp_bau['frac_org_fert_sd'].values[0])
        EF_vol = np.random.normal(temp_bau['ef_vol'].values[0], temp_bau['ef_vol_sd'].values[0])
        FRAC_LEACH = np.random.normal(temp_bau['frac_leach'].values[0], temp_bau['frac_leach_sd'].values[0])
        EF_leach = np.random.normal(temp_bau['ef_leaching'].values[0], temp_bau['ef_leaching_sd'].values[0])
        EF_PRP = np.random.normal(temp_bau['ef_prp'].values[0], temp_bau['ef_prp_sd'].values[0])
        dirN = (FSN + FON + FSOM) * EFdir + FPRP * EF_PRP  # Direct emissions
        volN = (FSN * FRAC_GASF + (FON + FPRP) * FRAC_GASM) * EF_vol  # Indirect volatilization
        leachN = (FSN + FON + FSOM + FPRP) * FRAC_LEACH * EF_leach  # indirect leaching
        N2O = (dirN + volN + leachN) / 1000 * 44 / 28  # sum and convert to tN2O

        results_unc[f"N2O_{aoi_id}"].append(N2O)
        
        #FIXME: add if incorporating fire 
        # Add change in N2O due to fire management
        # if "change in fire management" in intervention_subcategory:
        #     fire_n2o_ef = np.random.normal(temp_bau['burning_n2o_ef_mean'].values[0], temp_bau['burning_n2o_ef_sd'].values[0])
        #     # Add fires for bau scenario
        #     if temp_bau['fire_used'].values[0] == "True":
        #         CF_bau = np.random.normal(temp_bau['combustion_factor_mean'].values[0], temp_bau['combustion_factor_sd'].values[0])
        #         MB_bau = np.random.normal(temp_bau['fuel_biomass_mean'].values[0], temp_bau['fuel_biomass_sd'].values[0])
        #         fireN2O_bau = MB_bau * CF_bau * fire_n2o_ef / 1000
        #         N2O = N2O - fireN2O_bau
        #         fire_per_bau = temp_bau['fire_management_years'].values[0]
        #         fire_yrs_bau = [y for y in range(1, 20) if y % fire_per_bau == 0]
        #         results_unc[m, 4 + np.array(fire_yrs_bau)] = N2O - fireN2O_bau
        #     if temp_int['fire_used'].values[0] == "True":
        #         CF_int = np.random.normal(temp_int['combustion_factor_mean'].values[0], temp_int['combustion_factor_sd'].values[0])
        #         MB_int = np.random.normal(temp_int['fuel_biomass_mean'].values[0], temp_int['fuel_biomass_sd'].values[0])
        #         fireN2O_int = MB_int * CF_int * fire_n2o_ef / 1000
        #         results_unc[m, 4] += fireN2O_int
        #         fire_per_int = temp_int['fire_management_years'].values[0]
        #         fire_yrs_int = [y for y in range(1, 20) if y % fire_per_int == 0]
        #         results_unc[m, 4 + np.array(fire_yrs_int)] += fireN2O_int

        
        #FIXME: add if incorporating fire
        # # emissions from fire management
        # if "change in fire management" in intervention_subcategory:
        #     fire_ch4_ef = np.random.normal(temp_bau['burning_ch4_ef_mean'].values[0], temp_bau['burning_ch4_ef_sd'].values[0])
        #     # Add fires for bau scenario
        #     if temp_bau['fire_used'].values[0] == "True":
        #         fireCH4_bau = MB_bau * CF_bau * fire_ch4_ef / 1000
        #         results_unc[m, 24] = CH4_live - fireCH4_bau
        #         results_unc[m, 24 + np.array(fire_yrs_bau)] = CH4_live - fireCH4_bau
        #     if temp_int['fire_used'].values[0] == "True":
        #         fireCH4_int = MB_int * CF_int * fire_ch4_ef / 1000
        #         results_unc[m, 24] += fireCH4_int
        #         results_unc[m, 24 + np.array(fire_yrs_int)] += fireCH4_int


    # column_names = ['AOI', 'Area', 'Rep', 'SOC', 'totalC', f'N2O_y{aoi_id}', f'CH4_y{aoi_id}']
    results_unc_df = pd.DataFrame(results_unc)
    return results_unc_df


def generate_output():
    global biomass_co2_result
    json_data, df, AOIs, intervention_subcategory = load_json_data('FP_test.json')
    
    # Run Monte Carlo simulations
    np.random.seed(1)
    mc = [GHGcalc(aoi_id, df, 1000, intervention_subcategory, biomass_co2_result) for aoi_id in AOIs]


    # Process results
    result_list = []
    for iteration in mc:
        aoi_id = iteration['AOI'].iloc[0]  # Assuming 'AOI' is constant for each iteration
        result_dict = {'AOI': aoi_id}
        
        for column in iteration.columns:
            if column in ['AOI', 'Area', 'Rep']:
                result_dict[column] = iteration[column].iloc[0]  # Assuming these are constant for each iteration
            else:
                result_dict[column] = iteration[column].mean()
                result_dict[column+'_sd'] = iteration[column].std()
        
        result_list.append(result_dict)


    # Convert the list of dictionaries to a DataFrame
    result_df = pd.DataFrame(result_list)

    print(result_df.columns)

    # result_df has the following columns with aoi_id baked in
    # TODO: Need to roll them up for totals as well instead of aoi_id
    # Index(['AOI', 'Area', 'Rep', 'SOC', 'SOC_sd', 'totalC', 'totalC_sd',
    #    'N2O_tropicalmoist_acrisols', 'N2O_tropicalmoist_acrisols_sd',
    #    'CH4_tropicalmoist_acrisols', 'CH4_tropicalmoist_acrisols_sd',
    #    'N2O_tropicalmoist_ferralsols', 'N2O_tropicalmoist_ferralsols_sd',
    #    'CH4_tropicalmoist_ferralsols', 'CH4_tropicalmoist_ferralsols_sd'],
    #   dtype='object')
    
    # Summing specific columns
    result_df['N2O_thayr'] = result_df['N2O_tropicalmoist_acrisols'] + result_df['N2O_tropicalmoist_ferralsols']
    result_df['sdN2Oha'] = result_df['N2O_tropicalmoist_acrisols_sd'] + result_df['N2O_tropicalmoist_ferralsols_sd']
    result_df['CH4_thayr'] = result_df['CH4_tropicalmoist_acrisols'] + result_df['CH4_tropicalmoist_ferralsols']
    result_df['sdCH4ha'] = result_df['CH4_tropicalmoist_acrisols_sd'] + result_df['CH4_tropicalmoist_ferralsols_sd']
    result_df['CO2_thayr'] = result_df['SOC'] + result_df['totalC']
    result_df['sdCO2ha'] = result_df['SOC_sd'] + result_df['totalC_sd']
    result_df['CO2_tyr'] = result_df['CO2_thayr'] * result_df['Area']
    result_df['sdCO2'] = result_df['sdCO2ha'] * result_df['Area']
    result_df['N2O_tyr'] = result_df['N2O_thayr'] * result_df['Area']
    result_df['sdN2O'] = result_df['sdN2Oha'] * result_df['Area']
    result_df['CH4_tyr'] = result_df['CH4_thayr'] * result_df['Area']
    result_df['sdCH4'] = result_df['sdCH4ha'] * result_df['Area']

    # Drop the individual columns after summing
    result_df = result_df.drop(columns=[
        'SOC', 'totalC', 'SOC_sd', 'totalC_sd',
        'N2O_tropicalmoist_acrisols', 'N2O_tropicalmoist_ferralsols',
        'N2O_tropicalmoist_acrisols_sd', 'N2O_tropicalmoist_ferralsols_sd',
        'CH4_tropicalmoist_acrisols', 'CH4_tropicalmoist_ferralsols',
        'CH4_tropicalmoist_acrisols_sd', 'CH4_tropicalmoist_ferralsols_sd'
    ])

    # Aggregate the DataFrame by summing up the columns
    totals = result_df.sum(numeric_only=True).to_frame().transpose()
    totals['AOI'] = 'Total'
    totals['Rep'] = 'Total'
    result_df = result_df.append(totals, ignore_index=True)

    # Reorder columns to place 'AOI' first
    cols = ['AOI', 'Rep'] + [col for col in result_df.columns if col not in ['AOI', 'Rep']]
    result_df = result_df[cols]

    print(result_df)
    # 
    #NOTE: this is creating all the columns only with projections from year 1 up to you year 20
    column_names = ['AOI', 'Area', 'CO2_thayr'] + \
                   [f'N2O_thayr_y{i}' for i in range(1, 21)] + \
                   [f'CH4_thayr_y{i}' for i in range(1, 21)] + \
                   ['sdCO2ha'] + [f'sdN2Oha_y{i}' for i in range(1, 21)] + \
                   [f'sdCH4ha_y{i}' for i in range(1, 21)] + \
                   ['CO2_tyr'] + [f'N2O_tyr_y{i}' for i in range(1, 21)] + \
                   [f'CH4_tyr_y{i}' for i in range(1, 21)] + \
                   ['sdCO2'] + [f'sdN2O_y{i}' for i in range(1, 21)] + \
                   [f'sdCH4_y{i}' for i in range(1, 21)]

    #TODO: need to create a naming above
    resdf = pd.DataFrame(resdf, columns=column_names)

    #NOTE: setting the numbers for projections
    # # Convert to JSON output
    restib = []
    for r in range(len(resdf)):
        aoi_dict = {
            'AOI': f'Sub_AOI-{r+1}',
            'CO2_thayr': [resdf.iloc[r]['CO2_thayr']] * 20,
            # NOTE: biomass_co2_sd[aoi_id] gives the sd for biomass_co2_result for that aoi_id
            'sdCO2ha': [resdf.iloc[r]['sdCO2ha']] * 20,
            'CO2_tyr': [resdf.iloc[r]['CO2_tyr']] * 20,
            'sdCO2': [resdf.iloc[r]['sdCO2']] * 20,
            'N2O_thayr': resdf.iloc[r][[f'N2O_thayr_y{i}' for i in range(1, 21)]].tolist(),
            'sdN2Oha': resdf.iloc[r][[f'sdN2Oha_y{i}' for i in range(1, 21)]].tolist(),
            'N2O_tyr': resdf.iloc[r][[f'N2O_tyr_y{i}' for i in range(1, 21)]].tolist(),
            'sdN2O': resdf.iloc[r][[f'sdN2O_y{i}' for i in range(1, 21)]].tolist(),
            'CH4_thayr': resdf.iloc[r][[f'CH4_thayr_y{i}' for i in range(1, 21)]].tolist(),
            'sdCH4ha': resdf.iloc[r][[f'sdCH4ha_y{i}' for i in range(1, 21)]].tolist(),
            'CH4_tyr': resdf.iloc[r][[f'CH4_tyr_y{i}' for i in range(1, 21)]].tolist(),
            'sdCH4': resdf.iloc[r][[f'sdCH4_y{i}' for i in range(1, 21)]].tolist()
        }
        restib.append(aoi_dict)

    # Write the DataFrame to a CSV file
    result_df.to_csv('output\Forest_Protection_Deforestation_test_result.csv', index=False)
    # 
    # #TODO: need to create a json file from result_df 
    # json_output = json.dumps(result_df, indent=2)
    # 
    # #TODO: uncomment this part to create the file with santized output
    # # # Ensure the output directory exists
    # output_dir = 'output'
    # os.makedirs(output_dir, exist_ok=True)
    # 
    # # Write the JSON output
    # sanitized_filename = 'output_file.json'  # Replace with actual filename sanitization if needed
    # output_file = os.path.join(output_dir, sanitized_filename)
    # with open(output_file, 'w') as f:
    #     json.dump(json_data, f, indent=4)

    # sanitized = sanitize_filename(f'{intervention_subcategory}.json')
    # 
    # # # Write the JSON output
    # output_file = os.path.join(output_dir, sanitized)
    # with open(output_file, 'w') as f:
    #    f.write(json_output)

    #print(f"\n\nGenerated file: {sanitized_filename} in the output folder.")

generate_output()
