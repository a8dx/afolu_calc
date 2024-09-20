# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: Forest_Protection_Deforestation_test.py
# Author: Barbara Bomfim
# Date Started: 07/11/2024
# Last Edited: 09/20/2024
# Purpose: AFOLU GHG Calculations for Forest Protection from Deforestation (FP) Interventions
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data_FP from Dynamic World land use data_FP,  data_FP layers
#
# User inputs: User inputs depend on sub-intervention category (planted or natural forest; 
#  tillage; nutrient management; fire management). 
#   - Initial:
#           -land area: area under business-as-usual land use
#           -forest_type_deforestation: forest type prior to deforestation
#.          -forest_management_type: Natural, Plantation
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#           -fire_used: yes, no
#.  -Baseline data:
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
#  -Intervention requires:
#           -forest_area: forest area that is being protected
#           -forest_type: selection from dropdown list
#           -deforestation_rate_post: Deforestation rate (in %/yr) after intervention, based on published data_FP (e.g., global forest watch)
#           -fire_used: yes, no
#.          -fire_incidence_rate: default annual average fire incidence rate (0.1) but can be calculated as the average of the years of available data_FP. 
#               Calculate fraction of forest (by land cover data_FP) that had burned each year to get the annual fire incidence rate.
#           -bio_ef: gCO2 per kg of dry matter burnt

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
from typing import Dict, List, Any

import pandas as pd
import numpy as np
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor
import ee
import re
import math

from helpersFP import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock, sanitize_filename

### Step 0: troubleshooting
# Define function to check for required global variables
def check_global_variables():
    required_globals = ['data', 'polygon', 'average_agbd_dict']
    for var in required_globals:
        if var not in globals():
            raise NameError(f"Required global variable '{var}' is not defined")
          
check_global_variables()

### STEP 1 - Estimate mean annual AGB (in Mg/ha) using GEDI data_FPset ####
# GEE path: https://code.earthengine.google.com/?scriptPath=users%2Fbabomfimf%2FAGB-GEDI-L4A%3Atest-3_AGB_annual_mean
## Anthony's GEE: https://code.earthengine.google.com/44d6aaa5db4764b5b5f3825baf900d04

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
gedi_l4a_monthly = ee.ImageCollection('LARSE/GEDI/GEDI04_A_002_MONTHLY')

# Filter the collection to your polygon and a wider time range
filtered_gedi = gedi_l4a_monthly.filterBounds(polygon).filterDate('2019-01-01', '2020-01-01')#for bau, starts as far back as possible until just before intervention starts
# Filter agbd
agbd_band = filtered_gedi.select('agbd')
# Calculate annual mean
agbd_image = agbd_band.reduce(ee.Reducer.mean())

# Compute the mean annual AGBD value for the entire polygon
average_agbd = agbd_image.reduceRegion(
    reducer=ee.Reducer.mean(),
    geometry=polygon,
    scale=25  # Scale should match the resolution of the data_FPset
)
average_agbd_dict = average_agbd.getInfo()
print(f"{average_agbd_dict} Mg/ha")

average_agbd = agbd_image.reduceRegion(
    reducer=ee.Reducer.mean(),
    geometry=polygon,
    scale=25  # Scale should match the resolution of the data_FPset
)

print("\nArea of the polygon:")
print(f"{polygon.area().getInfo()/10000} ha") # to get area in hectares
maskedArea_ha = polygon.area().getInfo()/10000

## JSON input #####

# JSON file path - Get the current working directory
current_directory = os.getcwd()
print("Current working directory:", current_directory)
# Define the relative file path
file_name = 'FP_final.json'
file_path = os.path.join(current_directory, file_name)
json_file_path = './FP_final.json'

# Load data_FP from the JSON file
with open(file_path, 'r') as file:
        data = json.load(file)

### STEP 2: Calculate Mean annual Aboveground Carbon (AGC) stock in tCO2e/ha ####
# here we obtain aboveground biomass and convert to aboveground carbon stock
def mean_annual_agc_stock(scenario, log_level='info'): # edited this function
    global data, polygon, average_abgd_dict ## edited this line
    
    scenario_data = data["scenarios"][scenario][0]
    
    average_agbd_value = average_agbd_dict.get('agbd_mean', 0)
    area_converted = polygon.area().getInfo()/10000 # converting to hectares
    area_converted_yr = area_converted # area_converted_yr is defined
    
    CF = scenario_data.get("CF", 0) #dry matter to carbon conversion factor
    average_agbd_cstock = convert_to_c(average_agbd_value, CF)
    average_agbd_tco2e = convert_to_co2e(average_agbd_cstock)
    
    if log_level == 'debug':
        print(f"{area_converted} ha")
        print(f"Average Aboveground Biomass Density: {average_agbd_value} Mg/ha")
        print(f"Average Aboveground Carbon Stock: {average_agbd_tco2e} tCO2e/ha")
    
    subregion_results = []
    for subregion in scenario_data["aoi_subregions"]:
        subregion_area = subregion["area"]
        subregion_fraction = subregion_area / area_converted_yr if area_converted_yr !=0 else 0 # edited this line
        
        subregion_agbd_tco2e = average_agbd_tco2e * subregion_fraction
        
        subregion_results.append({
            "aoi_id": subregion["aoi_id"],
            "area": subregion_area,
            "carbon_stock": subregion_agbd_tco2e
        })
        
        if log_level == 'debug':
            print(f"Subregion {subregion['aoi_id']}: {subregion_agbd_tco2e} tCO2e")
            print(f"Mean annual AGC stock calculated for scenario: {scenario}")
            print(f"Total area converted: {area_converted_yr} ha")
            print(f"Average AGBD value: {average_agbd_value} Mg/ha")
            print(f"Average AGBD CO2e: {average_agbd_tco2e} tCO2e/ha")
    
    return subregion_results
  
### STEP 3 - Calculate Mean Annual Total C stock ####
# calculate the sum of aboveground and belowground carbon stocks
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

## Step 4.1: Define equation to calculate Historical Deforestation Rate 
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
  
## Step 4.2: Define equation to calculate Baseline Yearly Deforested Area
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
  
##Define equation to calculate Yearly Baseline Total Biomass Carbon Stock Loss
# def calculate_carbon_stock_loss_yeari(total_results, scenario, log_level='info'):
#     global data
#     
#     scenario_data = data["scenarios"][scenario][0]
#     #ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
#     
#     total_results = []
#     
#     for subregion in aoi_subregions:
#         aoi_id = subregion["aoi_id"]
#         deforested_area_yeari = subregion["deforested_area_year_i"]
#         average_total_tco2e = subregion["total_carbon_stock"]
#  
#         carbon_stock_loss = deforested_area_yeari * average_total_tco2e
#         
#         total_results.append({
#             "total_carbon_stock_loss": carbon_stock_loss
#         })
#         
#         if log_level == 'debug':
#             print(f"Subregion {total_results['aoi_id']}:")
#             print(f"  Carbon Stock Loss year i: {carbon_stock_loss:.2f} tCO2/yr")
#     
#     return total_results


#### Step 4.3: Define equation to calculate Baseline Annual Total Biomass CO2 Emissions ##

# def calculate_baseline_biomass_co2_emissions(total_results, scenario, log_level='info'):
#     global data
#     
#     scenario_data = data["scenarios"][scenario][0]
#     #ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
#     
#     total_results = []
#     
#     for subregion in aoi_subregions:
#         aoi_id = subregion["aoi_id"]
#         historical_def_rate = subregion["hist_def_rate"]
#         carbon_stock_loss_yeari = subregion["carbon_stock_loss"]
#  
#         annual_baseline_biomass_co2_emissions = historical_def_rate * carbon_stock_loss_yeari
#         
#         total_results.append({
#             "annual_baseline_biomass_co2": annual_baseline_biomass_co2_emissions
#         })
#         
#         if log_level == 'debug':
#             print(f"Subregion {total_results['aoi_id']}:")
#             print(f"  Annual Baseline Emissions year i: {annual_baseline_biomass_co2_emissions:.2f} tCO2/yr")
#     
#     return total_results

#edited code
def load_json_data(json_file_path):
    with open(json_file_path, 'r') as file:
        return json.load(file)

def calculate_historical_deforestation_rate(scenario_data):
    deforestation_pre = float(scenario_data.get("deforestation_rate_pre", 0))
    deforestation_post = float(scenario_data.get("deforestation_rate_post", 0))
    return deforestation_pre - deforestation_post

def calculate_deforested_area_yeari(forest_area, hist_def_rate):
    return (hist_def_rate / 100) * forest_area

def calculate_carbon_stock_loss_yeari(deforested_area, average_total_tco2e):
    return deforested_area * average_total_tco2e

def calculate_baseline_biomass_co2_emissions(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    
    years = int(scenario_data.get('years', 20))
    if years != 20:
        print(f"Warning: 'years' in the JSON file is set to {years}. Using 20 years for calculations.")
    years = 20  # Ensure we use 20 years for the calculation
    
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    hist_def_rate = calculate_historical_deforestation_rate(scenario_data)
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = float(subregion["area"])
        average_total_tco2e = float(subregion.get("total_carbon_stock", 0))
        
        subregion_results = {
            "aoi_id": aoi_id,
            "hist_def_rate": hist_def_rate,
            "annual_emissions": [],
            "total_emissions": 0,
            "total_deforested_area": 0
        }
        
        remaining_area = forest_area
        for year in range(1, years + 1):
            deforested_area_yeari = min(calculate_deforested_area_yeari(remaining_area, hist_def_rate), remaining_area)
            carbon_stock_loss_yeari = calculate_carbon_stock_loss_yeari(deforested_area_yeari, average_total_tco2e)
            annual_baseline_biomass_co2_emissions = carbon_stock_loss_yeari
            
            subregion_results["annual_emissions"].append({
                "year": year,
                "deforested_area": deforested_area_yeari,
                "carbon_stock_loss": carbon_stock_loss_yeari,
                "baseline_biomass_co2_emissions": annual_baseline_biomass_co2_emissions
            })
            
            subregion_results["total_emissions"] += annual_baseline_biomass_co2_emissions
            subregion_results["total_deforested_area"] += deforested_area_yeari
            
            remaining_area -= deforested_area_yeari
            
            if log_level == 'debug':
                print(f"Subregion {aoi_id}, Year {year}:")
                print(f"  Historical Deforestation Rate: {hist_def_rate:.2f}%/yr")
                print(f"  Deforested Area: {deforested_area_yeari:.2f} ha")
                print(f"  Carbon Stock Loss: {carbon_stock_loss_yeari:.2f} tCO2")
                print(f"  Annual Baseline Emissions: {annual_baseline_biomass_co2_emissions:.2f} tCO2/yr")
            
            if remaining_area <= 0:
                break
        
        total_results.append(subregion_results)
    
    return total_results

# Testing
json_file_path = 'FP_final.json'
scenario = 'business_as_usual'  # or 'intervention'

results = calculate_baseline_biomass_co2_emissions(json_file_path, scenario, log_level='debug')

if results:
    print("\nBaseline Biomass CO2 Emissions Results:")
    for subregion in results:
        print(f"AOI: {subregion['aoi_id']}")
        print(f"Historical Deforestation Rate: {subregion['hist_def_rate']:.2f}%/yr")
        print(f"  Total Deforested Area over 20 years: {subregion['total_deforested_area']:.2f} ha")
        print(f"  Total Baseline Biomass CO2 Emissions over 20 years: {subregion['total_emissions']:.2f} tCO2")
        
        if subregion['annual_emissions']:
            print("  Annual Emissions:")
            for year_data in subregion['annual_emissions']:
                print(f"    Year {year_data['year']}: {year_data['baseline_biomass_co2_emissions']:.2f} tCO2")
        else:
            print("  No annual emissions data available.")
else:
    print("Failed to calculate Baseline Biomass CO2 Emissions.")

# Print raw results for debugging
print("\nRaw Results:")
import pprint
pprint.pprint(results)




####### Intervention Calculations: Avoided Deforestation Total Biomass CO2 Calculations #######

### Step 5: Define equations needed to calculate Avoided Deforestation Emissions from Trees:
## Step 5.1: Equation to calculate Actual deforested area and Forest area end of year n
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

## Step 5.2: Calculate area of avoided deforestation in year n
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

## Step 5.3: Calculate Avoided Deforestation Emissions from Trees
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

### Step 6: Calculate Expected Carbon Stock changes due to Avoided Deforestation

##Step 6.1: Define function to estimate average annual biomass growth
def load_json_data(json_file_path):
    with open(json_file_path, 'r') as file:
        return json.load(file)

def calculate_average_annual_biomass_growth(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = float(subregion["area"])
        
        # Retrieve parameters from JSON data
        R = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))   # Root-to-shoot ratio
        Iv = float(scenario_data.get("Iv", 0))  # Average net annual increment (added to json file)
        BCEFi = float(scenario_data.get("BCEFi", 0))  # Biomass conversion and expansion factor (added to the json file)
        
        # Calculate Gtotal
        Gtotal = Iv * BCEFi * (1 + R)
        
        subregion_result = {
            "aoi_id": aoi_id,
            "area": area,
            "Gtotal": Gtotal,
            "parameters": {
                "R": R,
                "Iv": Iv,
                "BCEFi": BCEFi
            }
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion: {aoi_id}")
            print(f"  Area: {area:.2f} ha")
            print(f"  R: {R:.4f}")
            print(f"  Iv: {Iv:.4f} m3 ha-1 yr-1")
            print(f"  BCEFi: {BCEFi:.4f}")
            print(f"  Gtotal: {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
            print()
    
    return total_results

# Example usage
json_file_path = 'FP_final.json'
scenario = 'business_as_usual'  # or 'intervention'

results = calculate_average_annual_biomass_growth(json_file_path, scenario, log_level='debug')

if results:
    print("\nAverage Annual Biomass Growth Results:")
    for subregion in results:
        print(f"AOI: {subregion['aoi_id']}")
        print(f"  Area: {subregion['area']:.2f} ha")
        print(f"  Average Annual Biomass Growth (Gtotal): {subregion['Gtotal']:.4f} tonnes d.m. ha-1 yr-1")
        print("  Parameters used:")
        for param, value in subregion['parameters'].items():
            print(f"    {param}: {value:.4f}")
        print()
else:
    print("Failed to calculate Average Annual Biomass Growth.")

# Print raw results for debugging
print("\nRaw Results:")
import pprint
pprint.pprint(results)


##Step 6.2: Define function to estimate annual increase in biomass carbon stock due to biomass growth
def calculate_annual_biomass_carbon_stock_increase(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        A = float(subregion["area"])  # Area in hectares
        CF = float(scenario_data.get("CF", 0))  # Carbon fraction of dry matter
        
        # Calculate Gtotal using the previous function
        Gtotal = calculate_average_annual_biomass_growth(scenario_data, subregion)
        
        # Calculate ΔCG
        delta_CG = A * Gtotal * CF
        
        subregion_result = {
            "aoi_id": aoi_id,
            "area": A,
            "Gtotal": Gtotal,
            "CF": CF,
            "delta_CG": delta_CG
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion: {aoi_id}")
            print(f"  Area (A): {A:.2f} ha")
            print(f"  Average Annual Biomass Growth (Gtotal): {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
            print(f"  Carbon Fraction (CF): {CF:.4f} tonne C (tonne d.m.)-1")
            print(f"  Annual Increase in Biomass Carbon Stocks (ΔCG): {delta_CG:.4f} tC/yr")
            print()
    
    return total_results

# testing
json_file_path = 'FP_final.json'
scenario = 'business_as_usual'  # or 'intervention'

results = calculate_annual_biomass_carbon_stock_increase(json_file_path, scenario, log_level='debug')

if results:
    print("\nAnnual Increase in Biomass Carbon Stocks Results:")
    total_delta_CG = 0
    for subregion in results:
        print(f"AOI: {subregion['aoi_id']}")
        print(f"  Area: {subregion['area']:.2f} ha")
        print(f"  Average Annual Biomass Growth (Gtotal): {subregion['Gtotal']:.4f} tonnes d.m. ha-1 yr-1")
        print(f"  Carbon Fraction (CF): {subregion['CF']:.4f} tonne C (tonne d.m.)-1")
        print(f"  Annual Increase in Biomass Carbon Stocks (ΔCG): {subregion['delta_CG']:.4f} tC/yr")
        print()
        total_delta_CG += subregion['delta_CG']
    
    print(f"Total Annual Increase in Biomass Carbon Stocks: {total_delta_CG:.4f} tC/yr")
else:
    print("Failed to calculate Annual Increase in Biomass Carbon Stocks.")

# Print raw results for debugging
print("\nRaw Results:")
import pprint
pprint.pprint(results)

#### Step 6.3:Calculate Annual Decrease in Biomass Stock ####

## Estimate ∆CL using Equation 2.11: ∆CL = Lwood−removals + Lfuelwood + Ldisturbance
# Lwood-removals = annual carbon loss due to wood removals, tonnes C yr-1 (Equation 2.12) 
# Lfuelwood = annual biomass carbon loss due to fuelwood removals, tC/yr (Equation 2.13)
# Ldisturbance = annual biomass carbon losses due to disturbances, tonnes C yr-1 (See Equation 2.14)

## Step 6.3.1: Define Lwood-removals equation
def calculate_lwood_removals(scenario, log_level='debug'):
    """
    Calculate Lwood-removals using the provided equation for each subregion.
    Parameters:
    scenario (str): Scenario in question
    Returns:
    list: Annual carbon loss due to wood removals in tCO2e/yr for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    results = []
    total_lwood_removals = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        H = scenario_data.get("H")
        BCEFr = scenario_data.get("BCEFr") # biomass conversion and expansion factors applicable to wood removals; transforms
                                            # merchantable biomass to total biomass (including bark)
        ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        CF = scenario_data.get("CF")

        print(H, BCEFr, ratio, CF)
        
        if None in [H, BCEFr, ratio, CF]:
            if log_level == 'debug':
                print(f"Missing data for subregion {aoi_id} in scenario {scenario}. Please check your data.")
            continue
        
        lwood_removals = H * BCEFr * (1 + ratio) * CF * (44 / 12) * (area / scenario_data.get("area_converted_yr", 1))
        
        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "lwood_removals": lwood_removals
        })
        
        total_lwood_removals += lwood_removals
        total_area += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Annual carbon loss due to wood removals: {lwood_removals:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")
    
    average_lwood_removals = total_lwood_removals / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal annual carbon loss due to wood removals across all subregions: {total_lwood_removals:.2f} tCO2e/yr")
        print(f"Average annual carbon loss due to wood removals per hectare: {average_lwood_removals:.2f} tCO2e/ha/yr")
    
    return results

## Step 6.3.2:Define Lfuelwood equation
def calculate_lfuelwood(scenario, log_level='info'):
    """
    Calculate Lfuelwood using the provided equation for each subregion.
    Parameters:
    scenario (str): Scenario in question
    Returns:
    list: Annual biomass carbon loss due to fuelwood removals in tCO2e/yr for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    results = []
    total_lfuelwood = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        FGtrees = scenario_data.get("FGtrees")
        FGpart = scenario_data.get("FGpart")
        BCEFr = scenario_data.get("BCEFr")
        ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        D = scenario_data.get("D")
        CF = scenario_data.get("CF")
        
        if None in [FGtrees, FGpart, BCEFr, ratio, D, CF]:
            if log_level == 'debug':
                print(f"Missing data for calculating Lfuelwood in subregion {aoi_id}. Please check your data for scenario {scenario}.")
            continue
        
        lfuelwood = ((FGtrees * BCEFr * (1 + ratio)) + FGpart * D) * CF * (44 / 12) * (area / scenario_data.get("area_converted_yr", 1))
        
        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "lfuelwood": lfuelwood
        })
        
        total_lfuelwood += lfuelwood
        total_area += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Annual biomass carbon loss due to fuelwood removals: {lfuelwood:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")
    
    average_lfuelwood = total_lfuelwood / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal annual biomass carbon loss due to fuelwood removals across all subregions: {total_lfuelwood:.2f} tCO2e/yr")
        print(f"Average annual biomass carbon loss due to fuelwood removals per hectare: {average_lfuelwood:.2f} tCO2e/ha/yr")
    
    return results

## Step 6.3.3: Estimate Annual Decrease in Biomass Carbon Stock
#Define Ldisturbance equation
def calculate_ldisturbance(scenario, log_level='info'):
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    results = []
    total_ldisturbance = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        Adisturbance = scenario_data.get("Adisturbance")
        Bw = scenario_data.get("Bw")
        ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        fd = scenario_data.get("fd")
        
        if None in [Adisturbance, Bw, ratio, fd]:
            if log_level == 'debug':
                print(f"Missing data for calculating Ldisturbance in subregion {aoi_id}. Please check your data.")
            continue
        
        ldisturbance = Adisturbance * Bw * (1 + ratio) * fd * (44 / 12) * (area / scenario_data.get("area_converted_yr", 1))
        
        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "ldisturbance": ldisturbance
        })
        
        total_ldisturbance += ldisturbance
        total_area += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Annual biomass carbon losses due to disturbances: {ldisturbance:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")
    
    average_ldisturbance = total_ldisturbance / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal annual biomass carbon losses due to disturbances across all subregions: {total_ldisturbance:.2f} tCO2e/yr")
        print(f"Average annual biomass carbon losses due to disturbances per hectare: {average_ldisturbance:.2f} tCO2e/ha/yr")
    
    return results

## Step 6.3.4:Define ∆CL equation (∆CL = Lwood −removals + Lfuelwood + Ldisturbance)
def calculate_delta_cl(lwood_removals_results, lfuelwood_results, ldisturbance_results, log_level='info'):
    """
    Calculate the change in biomass stock (∆CL) for each subregion.
    Parameters:
    lwood_removals_results (list): List of dictionaries containing Lwood-removals results for each subregion.
    lfuelwood_results (list): List of dictionaries containing Lfuelwood results for each subregion.
    ldisturbance_results (list): List of dictionaries containing Ldisturbance results for each subregion.
    log_level (str): Logging level. Set to 'debug' for detailed output.
    Returns:
    list: Change in biomass stock (∆CL) in tCO2e/yr for each subregion.
    """
    results = []
    total_delta_cl = 0
    total_area = 0


    print("<<<<<>>>>", lwood_removals_results, lfuelwood_results, ldisturbance_results)

    # Ensure all input lists have the same length
    if not (len(lwood_removals_results) == len(lfuelwood_results) == len(ldisturbance_results)):
        if log_level == 'debug':
            print("Error: Input lists have different lengths. Please ensure all subregions are represented in each input.")
        return None

    for lwood, lfuel, ldist in zip(lwood_removals_results, lfuelwood_results, ldisturbance_results):
        # Ensure we're dealing with the same subregion in all inputs
        if not (lwood['aoi_id'] == lfuel['aoi_id'] == ldist['aoi_id']):
            if log_level == 'debug':
                print(f"Error: Mismatched subregion IDs: {lwood['aoi_id']}, {lfuel['aoi_id']}, {ldist['aoi_id']}")
            continue

        aoi_id = lwood['aoi_id']
        area = lwood['area']  # Assuming area is the same in all inputs

        # Calculate ∆CL for the subregion
        delta_cl = lwood['lwood_removals'] + lfuel['lfuelwood'] + ldist['ldisturbance']

        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "delta_cl": delta_cl
        })

        total_delta_cl += delta_cl
        total_area += area

        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Change in biomass stock (∆CL): {delta_cl:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")

    average_delta_cl = total_delta_cl / total_area if total_area > 0 else 0

    if log_level == 'debug':
        print(f"\nTotal change in biomass stock (∆CL) across all subregions: {total_delta_cl:.2f} tCO2e/yr")
        print(f"Average change in biomass stock (∆CL) per hectare: {average_delta_cl:.2f} tCO2e/ha/yr")

    print("RESSSSS", results) # note to edit this line or comment it

    return results  

### Step 6.4: Estimate avoided biomass burning emissions from avoided fire
def calculate_avoided_biomass_burning_emissions(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    Cf = 0.34  # Combustion factor, Value obtained from winrock AFOLU tool
    bio_ef = float(scenario_data.get("bio_ef", 0)) #from json file, (1580 g CO2 per kg of dry matter burnt)
    
    # Retrieve fire incidence rate from JSON
    fire_incidence_rate = float(scenario_data.get("fire_incidence_rate", 0))
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = float(subregion["area"])  # Forest area in hectares
        tree_carbon_stock = float(subregion.get("tree_carbon_stock", 0))  # Tree Carbon Stock in tC/ha
        
        # Calculate avoided biomass burning emissions
        avoided_emissions = (
            forest_area * 
            (fire_incidence_rate * tree_carbon_stock * (1/0.47) * Cf * bio_ef) * 
            1e-3  # Convert to tCO2/yr
        )
        
        subregion_result = {
            "aoi_id": aoi_id,
            "forest_area": forest_area,
            "tree_carbon_stock": tree_carbon_stock,
            "avoided_biomass_burning_emissions": avoided_emissions
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion: {aoi_id}")
            print(f"  Forest Area: {forest_area:.2f} ha")
            print(f"  Tree Carbon Stock: {tree_carbon_stock:.2f} tC/ha")
            print(f"  Fire Incidence Rate: {fire_incidence_rate:.4f}")
            print(f"  Combustion Factor (Cf): {Cf:.2f}")
            print(f"  Fire Emissions Factor (bio_ef): {bio_ef:.4f}")
            print(f"  Avoided Biomass Burning Emissions: {avoided_emissions:.4f} tCO2/yr")
            print()
    
    return total_results

# Example usage
json_file_path = 'IFM_final.json'
scenario = 'business_as_usual'  # or 'intervention'

results = calculate_avoided_biomass_burning_emissions(json_file_path, scenario, log_level='debug')

if results:
    print("\nAvoided Fire and Biomass Burning Emissions Results:")
    total_avoided_emissions = 0
    for subregion in results:
        print(f"AOI: {subregion['aoi_id']}")
        print(f"  Forest Area: {subregion['forest_area']:.2f} ha")
        print(f"  Tree Carbon Stock: {subregion['tree_carbon_stock']:.2f} tC/ha")
        print(f"  Avoided Biomass Burning Emissions: {subregion['avoided_biomass_burning_emissions']:.4f} tCO2/yr")
        print()
        total_avoided_emissions += subregion['avoided_biomass_burning_emissions']
    
    print(f"Total Avoided Biomass Burning Emissions: {total_avoided_emissions:.4f} tCO2/yr")
else:
    print("Failed to calculate Avoided Fire and Biomass Burning Emissions.")

# Print raw results for debugging
print("\nRaw Results:")
import pprint
pprint.pprint(results)


### TO DISCUSS ###


#### Step 6.5: Estimate Emissions from Illegal Logging

### Step 6.5.1: Emissions from Incidental Damage from Illegal Logging

######## question: how do we obtain illegal logging rate (ILR)??? Need to include in json input file.



##Define equation to calculate Carbon Emissions from Incidental Damage
def calculate_incidental_damage_illegal(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data["aoi_subregions"]
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        VolExt = float(scenario_data.get("VolExt", 0))
        LDF_factor_1 = float(scenario_data.get("LDF_factor_1", 0))
        LDF_factor_2 = float(scenario_data.get("LDF_factor_2", 0))
        ILR = float(scenario_data.get("ILR", 0))
        average_total_tco2e = float(scenario_data.get("average_total_tco2e", 0))
        
        LDF = LDF_factor_1 * average_total_tco2e + LDF_factor_2
        IncidentalDamage_illegal = LDF * ILR
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "LDF": LDF,
            "IncidentalDamage_illegal": IncidentalDamage_illegal
        })
    
    return total_results

##Define equation to calculate Illegal Logging Emissions from Trees
def calculate_logging_emissions_illegal(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        #AHA = float(scenario_data.get("AHA", 0)) # this should be user input, as mentioned in User Input commented part of python code
        #VolExt = float(scenario_data.get("VolExt", 0))
        D = float(scenario_data.get("D", 0))
        ElE_factor_1 = float(scenario_data.get("ElE_factor_1", 0))
        ElE_factor_2 = float(scenario_data.get("ElE_factor_2", 0))
        
        ELE = (ElE_factor_1 * D) - ElE_factor_2
        TimberTree_illegal = area * ILR * ELE
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "TimberTree": TimberTree,
            "ELE": ELE
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Extracted Log Emissions (ELE): {ELE:.4f} tC/m3 extracted")
            print(f"  Illegal Logging emissions : {TimberTree_illegal:.2f} tC")
    
    return total_results

##### Step 6.5.2: Estimate emissions from Community offtake
def calculate_community_offtake_emissions(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        community_area = float(subregion.get("area", 0))
        community_offtake = float(subregion.get("community_offtake", 0))
        
        D = float(scenario_data.get("D", 0))  # Wood Density
        ElE_factor_1 = float(scenario_data.get("ElE_factor_1", 0))
        ElE_factor_2 = float(scenario_data.get("ElE_factor_2", 0))
        
        # Calculate ELE
        ELE = (ElE_factor_1 * D) - ElE_factor_2
        
        # Calculate TimberTreecommunity
        TimberTreecommunity = community_area * community_offtake * ELE
        
        # Calculate Incidental DamageComm (assumed to be a fraction of TimberTreecommunity)
        # You may need to adjust this calculation based on your specific requirements
        incidental_damage_fraction = float(scenario_data.get("incidental_damage_fraction", 0.1))
        IncidentalDamageComm = TimberTreecommunity * incidental_damage_fraction
        
        # Calculate Community Offtake emissions
        community_offtake_emissions = (TimberTreecommunity + IncidentalDamageComm) * 44/12
        
        subregion_result = {
            "aoi_id": aoi_id,
            "community_area": community_area,
            "community_offtake": community_offtake,
            "ELE": ELE,
            "TimberTreecommunity": TimberTreecommunity,
            "IncidentalDamageComm": IncidentalDamageComm,
            "community_offtake_emissions": community_offtake_emissions
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Community Area: {community_area:.2f} ha")
            print(f"  Community Offtake: {community_offtake:.4f}")
            print(f"  Extracted Log Emissions (ELE): {ELE:.4f} tC/m3 extracted")
            print(f"  TimberTreecommunity: {TimberTreecommunity:.2f} tC")
            print(f"  Incidental DamageComm: {IncidentalDamageComm:.2f} tC")
            print(f"  Community Offtake Emissions: {community_offtake_emissions:.2f} tCO2/yr")
            print()
    
    return total_results

# Testing
json_file_path = 'FP_final.json'
scenario = 'business_as_usual'  # or 'intervention'

results = calculate_community_offtake_emissions(json_file_path, scenario, log_level='debug')

if results:
    print("\nCommunity Offtake Emissions Results:")
    total_emissions = 0
    for subregion in results:
        print(f"AOI: {subregion['aoi_id']}")
        print(f"  Community Area: {subregion['community_area']:.2f} ha")
        print(f"  Community Offtake: {subregion['community_offtake']:.4f}")
        print(f"  Community Offtake Emissions: {subregion['community_offtake_emissions']:.2f} tCO2/yr")
        print()
        total_emissions += subregion['community_offtake_emissions']
    
    print(f"Total Community Offtake Emissions: {total_emissions:.2f} tCO2/yr")
else:
    print("Failed to calculate Community Offtake Emissions.")

# Print raw results for debugging
print("\nRaw Results:")
import pprint
pprint.pprint(results)


####### Step 7: Calculate all - equations for baseline and intervention emissions calculation

### Step 7.1: Calculate all - baseline emissions
def calculate_all_bau(scenario, log_level='info'):
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
result["yearly_baseline_emissions_trees"]


### Step 7.2: Calculate all - equations for intervention emissions calculation
def calculate_all_bau(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    actual_area_deforested_each_year = calculate_actual_area_deforested_year_n(scenario=scenario, log_level=log_level)
    result["actual_area_deforested_each_year"] = actual_area_deforested_each_year
    
    forest_area_each_year = calculate_forest_area_end_of_year_n(scenario=scenario, n=20)
    result["forest_area_each_year"] = forest_area_each_year
    
    area_avoided_deforestation_year_n = calculate_area_avoided_deforestation_year_n(scenario=scenario, log_level=log_level)
    result["area_avoided_deforestation_year_n"] = area_avoided_deforestation_year_n
    
    yearly_avoided_emissions_trees = calculate_avoided_emissions_trees_each_year(scenario=scenario, log_level=log_level)
    result["yearly_avoided_emissions_trees"] = yearly_avoided_emissions_trees
    
    yearly_biomass_growth = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["yearly_biomass_growth"] = yearly_biomass_growth
    
    yearly_increase_biomass_stock = calculate_annual_biomass_carbon_stock_increase(scenario=scenario, log_level=log_level)
    result["yearly_increase_biomass_stock"] = yearly_increase_biomass_stock
    
    yearly_decrease_biomass_stock= calculate_delta_cl(scenario=scenario, log_level=log_level)
    result["yearly_decrease_biomass_stock"] = yearly_decrease_biomass_stock
    
    yearly_avoided_burning = calculate_avoided_biomass_burning_emissions(scenario=scenario, log_level=log_level)
    result["yearly_decrease_biomass_stock"] = yearly_decrease_biomass_stock
    
    loss_incidental_damage_illegal = calculate_incidental_damage_illegal(scenario=scenario, log_level=log_level)
    result["incidental_damage_illegal"] = incidental_damage_illegal
    
    loss_illegal_logging_emissions = calculate_logging_emissions_illegal(scenario=scenario, log_level=log_level)
    result["loss_illegal_logging_emissions"] = loss_illegal_logging_emissions
    
    community_offtake_emissions = calculate_community_offtake_emissions(scenario=scenario, log_level=log_level)
    result["community_offtake_emissions"] = community_offtake_emissions
    
    return result
result["yearly_avoided_emissions_trees"]

#### Step 8: Calculate bau and intervention scenarios and difference between both

biomass = {
    "business_as_usual": calculate_all_bau(scenario="business_as_usual"),
    "intervention": calculate_all_int(scenario="intervention")
}

biomass_co2_result = {}
biomass_co2_sd = {}
biomass_co2_result_error_positive = 10  # this is a percentage

def get_annual_change(scenario_result):
    if 'annual_change_in_carbon_stocks' in scenario_result:
        return scenario_result['annual_change_in_carbon_stocks']
    else:
        print(f"Warning: 'annual_change_in_carbon_stocks' not found in {scenario_result.get('scenario', 'unknown')} scenario.")
        return []

intervention_changes = get_annual_change(biomass["intervention"])
bau_changes = get_annual_change(biomass["business_as_usual"])

for inter_item, bau_item in zip(intervention_changes, bau_changes):
    if not (isinstance(inter_item, dict) and isinstance(bau_item, dict)):
        print(f"Warning: Invalid item format. Skipping. Inter: {inter_item}, BAU: {bau_item}")
        continue

    if 'aoi_id' not in inter_item or 'aoi_id' not in bau_item:
        print(f"Warning: Missing aoi_id. Skipping. Inter: {inter_item}, BAU: {bau_item}")
        continue

    if inter_item['aoi_id'] != bau_item['aoi_id']:
        raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
    
    aoi_id = inter_item['aoi_id']
    
    if 'annual_change' not in inter_item or 'annual_change' not in bau_item:
        print(f"Warning: Missing annual_change for aoi_id {aoi_id}. Skipping.")
        continue

    difference = inter_item['annual_change'] - bau_item['annual_change']
    
    biomass_co2_result[aoi_id] = difference
    biomass_co2_sd[aoi_id] = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96

print("Biomass CO2 Result:", biomass_co2_result)
print("Biomass CO2 Standard Deviation:", biomass_co2_sd)


######STEP 9: SOIL GHG CALCULATIONS #########

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

#########################################################

## STEP 10
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
    # NOTE: these are not required with the aoi_id proposed
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

## Step 11: Define function to organize results

def sanitize_filename(filename: str) -> str:
    # Remove invalid characters
    sanitized = re.sub(r'[<>:"/\\|?*]', '', filename)
    # Replace spaces with underscores
    sanitized = sanitized.replace(' ', '_')
    # Limit length to 255 characters
    return sanitized[:255]

def average_results_sims(input_dict: Dict[str, Any]) -> Dict[str, Any]:
    results_unc = {}
    for key, value in input_dict.items():
        if isinstance(value, list):
            results_unc[key] = np.mean(value) if value else 0
        else:
            results_unc[key] = value
    return results_unc

### Step 10.3: Define ghg function
def GHGcalc(aoi_id: str, df: pd.DataFrame, nx: int, intervention_subcategory: str, biomass_co2_result: Dict[str, float]) -> pd.DataFrame:
    temp_bau = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "business-as-usual")]
    temp_int = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "intervention")]
    
    if temp_bau.empty or temp_int.empty:
        print(f"Warning: No data found for AOI {aoi_id}")
        return pd.DataFrame()
    
    results_unc = {
        'AOI': aoi_id,
        'Area': temp_bau['area'].values[0],
        'SOC': [],
        'totalC': [],
        f'N2O_{aoi_id}': [],
        f'CH4_{aoi_id}': []
    }
    
    for _ in range(nx):
        # SOC Stock Change calculation
        SOCREF = np.random.normal(temp_bau['SOC_ref_tonnes_C_ha'].values[0], temp_bau['SOC_ref_tonnes_C_ha_sd'].values[0])
        FLUbau, FMGbau, FIbau = temp_bau[['FLU', 'FMG', 'FI']].values[0]
        FLUint, FMGint, FIint = temp_int[['FLU', 'FMG', 'FI']].values[0]
        
        SOCbau = SOCREF * FLUbau * np.random.normal(FMGbau, temp_bau['FMG_sd'].values[0]) * np.random.normal(FIbau, temp_bau['FI_sd'].values[0])
        SOCint = SOCREF * FLUint * np.random.normal(FMGint, temp_int['FMG_sd'].values[0]) * np.random.normal(FIint, temp_int['FI_sd'].values[0])
        dSOC = SOCbau - SOCint
        results_unc["SOC"].append(dSOC)
        
        # Total carbon calculation
        biomass_co2 = biomass_co2_result.get(aoi_id, 0)  # Use 0 if aoi_id is not in the dictionary
        results_unc["totalC"].append(dSOC * 44/12 + biomass_co2)

        # N2O emissions calculation
        if "nutrient management" in intervention_subcategory:
            FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
                   temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100)
            FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
                   temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100)
        else:
            FSN = FON = 0
        
        FSOM = dSOC / 10 * 1000
        
        EFdir = np.random.normal(temp_bau['ef_direct_n2o'].values[0], temp_bau['ef_direct_n2o_sd'].values[0])
        FRAC_GASF = np.random.normal(temp_bau['frac_syn_fert'].values[0], temp_bau['frac_syn_fert_sd'].values[0])
        FRAC_GASM = np.random.normal(temp_bau['frac_org_fert'].values[0], temp_bau['frac_org_fert_sd'].values[0])
        EF_vol = np.random.normal(temp_bau['ef_vol'].values[0], temp_bau['ef_vol_sd'].values[0])
        FRAC_LEACH = np.random.normal(temp_bau['frac_leach'].values[0], temp_bau['frac_leach_sd'].values[0])
        EF_leach = np.random.normal(temp_bau['ef_leaching'].values[0], temp_bau['ef_leaching_sd'].values[0])
        
        dirN = (FSN + FON + FSOM) * EFdir
        volN = (FSN * FRAC_GASF + FON * FRAC_GASM) * EF_vol
        leachN = (FSN + FON + FSOM) * FRAC_LEACH * EF_leach
        N2O = (dirN + volN + leachN) / 1000 * 44 / 28
        
        results_unc[f"N2O_{aoi_id}"].append(N2O)

        # CH4 emissions calculation
        CH4_live = 0  # Initialize CH4 emissions

        # Fire management calculations
        if "change in fire management" in intervention_subcategory:
            fire_n2o_ef = np.random.normal(temp_bau['burning_n2o_ef_mean'].values[0], temp_bau['burning_n2o_ef_sd'].values[0])
            fire_ch4_ef = np.random.normal(temp_bau['burning_ch4_ef_mean'].values[0], temp_bau['burning_ch4_ef_sd'].values[0])
            #co equation
            #nox equation
            
            # BAU scenario
            if temp_bau['fire_used'].values[0] == "True":
                CF_bau = np.random.normal(temp_bau['combustion_factor_mean'].values[0], temp_bau['combustion_factor_sd'].values[0])
                MB_bau = np.random.normal(temp_bau['fuel_biomass_mean'].values[0], temp_bau['fuel_biomass_sd'].values[0])
                fireN2O_bau = MB_bau * CF_bau * fire_n2o_ef / 1000
                fireCH4_bau = MB_bau * CF_bau * fire_ch4_ef / 1000
                N2O -= fireN2O_bau
                CH4_live -= fireCH4_bau
            
            # Intervention scenario
            if temp_int['fire_used'].values[0] == "True":
                CF_int = np.random.normal(temp_int['combustion_factor_mean'].values[0], temp_int['combustion_factor_sd'].values[0])
                MB_int = np.random.normal(temp_int['fuel_biomass_mean'].values[0], temp_int['fuel_biomass_sd'].values[0])
                fireN2O_int = MB_int * CF_int * fire_n2o_ef / 1000
                fireCH4_int = MB_int * CF_int * fire_ch4_ef / 1000
                N2O += fireN2O_int
                CH4_live += fireCH4_int

        results_unc[f"CH4_{aoi_id}"].append(CH4_live)

    return pd.DataFrame(results_unc)
  
### Step 10.4: Define function to generate output

## Define equation to generate output
def generate_output(input_json: str):
    global biomass_co2_result
    json_data, df, AOIs, intervention_subcategory = load_json_data(input_json)

    np.random.seed(1)
    mc = [GHGcalc(aoi_id, df, 1000, intervention_subcategory, biomass_co2_result) for aoi_id in AOIs]
    mc = [result for result in mc if not result.empty]  # to filter out empty dataframes

    if not mc:
        print("Error: No valid results generated. Check input data.")
        return None

    result_list = []
    for iteration in mc:
        aoi_id = iteration['AOI'].iloc[0]
        result_dict = {'AOI': aoi_id, 'Area': iteration['Area'].iloc[0]}

        for column in iteration.columns:
            if column not in ['AOI', 'Area']:
                result_dict[column] = iteration[column].mean()
                result_dict[f"{column}_sd"] = iteration[column].std()

        result_list.append(result_dict)

    result_df = pd.DataFrame(result_list)

    # Calculate totals
    n2o_columns = [col for col in result_df.columns if col.startswith('N2O_') and not col.endswith('_sd')]
    ch4_columns = [col for col in result_df.columns if col.startswith('CH4_') and not col.endswith('_sd')]
    
    if n2o_columns:
        result_df['N2O_thayr'] = result_df[n2o_columns].sum(axis=1)
        result_df['sdN2Oha'] = np.sqrt((result_df[[f"{col}_sd" for col in n2o_columns]]**2).sum(axis=1))
    else:
        result_df['N2O_thayr'] = 0
        result_df['sdN2Oha'] = 0
    
    if ch4_columns:
        result_df['CH4_thayr'] = result_df[ch4_columns].sum(axis=1)
        result_df['sdCH4ha'] = np.sqrt((result_df[[f"{col}_sd" for col in ch4_columns]]**2).sum(axis=1))
    else:
        result_df['CH4_thayr'] = 0
        result_df['sdCH4ha'] = 0
    
    if 'SOC' in result_df.columns and 'totalC' in result_df.columns:
        result_df['CO2_thayr'] = result_df['SOC'] + result_df['totalC']
        result_df['sdCO2ha'] = np.sqrt(result_df['SOC_sd']**2 + result_df['totalC_sd']**2)
    else:
        result_df['CO2_thayr'] = 0
        result_df['sdCO2ha'] = 0

    result_df['CO2_tyr'] = result_df['CO2_thayr'] * result_df['Area']
    result_df['sdCO2'] = result_df['sdCO2ha'] * result_df['Area']
    result_df['N2O_tyr'] = result_df['N2O_thayr'] * result_df['Area']
    result_df['sdN2O'] = result_df['sdN2Oha'] * result_df['Area']
    result_df['CH4_tyr'] = result_df['CH4_thayr'] * result_df['Area']
    result_df['sdCH4'] = result_df['sdCH4ha'] * result_df['Area']

    # Drop individual columns after totaling
    columns_to_keep = ['AOI', 'Area', 'CO2_thayr', 'sdCO2ha', 'CO2_tyr', 'sdCO2', 
                       'N2O_thayr', 'sdN2Oha', 'N2O_tyr', 'sdN2O', 
                       'CH4_thayr', 'sdCH4ha', 'CH4_tyr', 'sdCH4']
    result_df = result_df[columns_to_keep]

    # Add total row
    totals = result_df.sum(numeric_only=True).to_frame().T
    totals['AOI'] = 'Total'
    result_df = pd.concat([result_df, totals], ignore_index=True)

    # Prepare JSON output
    json_output = []
    for _, row in result_df.iterrows():
        aoi_dict = {
            'AOI': row['AOI'],
            'Area': float(row['Area']),
            'CO2_thayr': [float(row['CO2_thayr'])] * 20,
            'sdCO2ha': [float(row['sdCO2ha'])] * 20,
            'CO2_tyr': [float(row['CO2_tyr'])] * 20,
            'sdCO2': [float(row['sdCO2'])] * 20,
            'N2O_thayr': [float(row['N2O_thayr'])] * 20,
            'sdN2Oha': [float(row['sdN2Oha'])] * 20,
            'N2O_tyr': [float(row['N2O_tyr'])] * 20,
            'sdN2O': [float(row['sdN2O'])] * 20,
            'CH4_thayr': [float(row['CH4_thayr'])] * 20,
            'sdCH4ha': [float(row['sdCH4ha'])] * 20,
            'CH4_tyr': [float(row['CH4_tyr'])] * 20,
            'sdCH4': [float(row['sdCH4'])] * 20
        }
        json_output.append(aoi_dict)

    # Generate output files
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)

    # Write CSV output (optional)
    csv_file = os.path.join(output_dir, 'result.csv')
    result_df.to_csv(csv_file, index=False)
    print(f"Generated CSV file: {csv_file}")

    # Write JSON output
    sanitized_filename = sanitize_filename(f'{intervention_subcategory}.json')
    json_file = os.path.join(output_dir, sanitized_filename)
    with open(json_file, 'w') as f:
        json.dump(json_output, f, indent=2)

    print(f"Generated JSON file: {json_file}")

    # Test if the JSON file was generated correctly
    if os.path.exists(json_file):
        with open(json_file, 'r') as f:
            try:
                test_load = json.load(f)
                print(f"JSON file successfully generated and can be loaded.")
            except json.JSONDecodeError:
                print(f"Error: JSON file was created but contains invalid JSON.")
    else:
        print(f"Error: JSON file was not created at {json_file}")

    return json_file

# Testing
try:
    json_file_path = './ARR_final.json'
    output_file = generate_output(json_file_path)
    if output_file:
        print(f"Output file saved at: {output_file}")
    else:
        print("Failed to generate output file.")
except Exception as e:
    print(f"An error occurred: {str(e)}")

## end ###
