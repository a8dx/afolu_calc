# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: Forest_Protection_Deforestation_test.py
# Author: Barbara Bomfim
# Date Started: 07/11/2024
# Last Edited: 09/29/2024
# Purpose: AFOLU GHG Calculations for Forest Protection from Deforestation (FP) Interventions
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data_FP from Dynamic World land use data_FP,  data_FP layers
#
# User inputs: User inputs depend on sub-intervention category: 
#               protection_from_fire, 
#               protection_from_deforestation,
#               protection_from_illegal_logging. 
#   - Initial:
#           -land area: area under business-as-usual land use
#           -forest_type_deforestation: forest type prior to deforestation
#.          -forest_management_type: Natural, Plantation
#           -fire_used: yes, no
#.  -Baseline data:
#           -deforested_land_area: area under business-as-usual land use (deforestation)
#           -deforestation_rate_pre: Deforestation rate (in %/yr) before intervention, based on published data_FP (e.g., global forest watch)
#           -forest_type_deforestation: forest type prior to deforestation
#.          -forest_management_type: Natural, Plantation
#           -fire_used: yes, no
#           -fire_frequency = 5 years, 10 years
#           -Bw = average AGB of land areas affected by disturbance (t dm/ha) is the same as the AGB from the GEDI. 
#           -fd = fraction of biomass lost in disturbance (dimensionless; default 0.34 for fire))
#           -illegal_logging_rate: volume timber over bark extracted illegally each year (m3 ha-1 yr-1)
#           -Gw: average annual increment in AGB in natural regeneration (t dm/ha/yr) - IPCC Tables 3A.1.5 and 3A.1.6
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
def mean_annual_agc_stock(scenario, log_level='info'):
    global data, average_agbd_dict
    
    scenario_data = data["scenarios"][scenario][0]
    
    average_agbd_value = average_agbd_dict.get('agbd_mean', 0)
    
    CF = scenario_data.get("CF", 0)
    average_agbd_cstock = convert_to_c(average_agbd_value, CF)
    average_agbd_tco2e = convert_to_co2e(average_agbd_cstock)
    
    subregion_results = []
    for subregion in scenario_data["aoi_subregions"]:
        forest_area = subregion["area"]
        
        subregion_agbd_tco2e = average_agbd_tco2e * forest_area
        
        subregion_results.append({
            "aoi_id": subregion["aoi_id"],
            "area": forest_area,
            "carbon_stock": subregion_agbd_tco2e
        })
        
        if log_level == 'debug':
            print(f"Subregion {subregion['aoi_id']}:")
            print(f"  Forest Area: {forest_area:.2f} ha")
            print(f"  Average AGBD: {average_agbd_value:.2f} Mg/ha")
            print(f"  Carbon Fraction (CF): {CF:.2f}")
            print(f"  Average AGBD Carbon Stock: {average_agbd_cstock:.2f} tC/ha")
            print(f"  Average AGBD CO2e: {average_agbd_tco2e:.2f} tCO2e/ha")
            print(f"  Subregion Carbon Stock: {subregion_agbd_tco2e:.2f} tCO2e")
            print()
    
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


####STEP 4A - PROTECTION FROM DEFORESTATION SUB-INTERVENTION ########
### BASELINE: Protection from Deforestation Annual Total Biomass CO2 calculations ###

## Equations:
# historical_deforestation_rate = calculate_historical_deforestation_rate(deforestation_rate_pre, deforestation_rate_post)
# area_deforested_value = calculate_area_deforested(historical_deforestation_rate, forest_area)
# carbon_stock_loss = calculate_carbon_stock_loss(area_deforested_value, average_total_tco2e)
# annual_baseline_biomass_co2_emissions = calculate_annual_baseline_emissions(historical_deforestation_rate, carbon_stock_loss)

## Step 4A.1: Define equation to calculate Historical Deforestation Rate 
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
  
## Step 4A.2: Define equation to calculate Baseline Yearly Deforested Area
def calculate_deforested_area_yeari(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])

    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area_aoi = subregion["area"] ## Anthony here check for forest_area (the area that is actually forest)
        historical_def_rate_aoi = subregion["hist_def_rate"]
    
    deforested_area_yeari = (historical_def_rate_aoi / 100) * forest_area_aoi
    
    total_results.append({
        "deforested_area_year_i": deforested_area_yeari
    })
    
    if log_level == 'debug':
        # print(f"Subregion {total_results['aoi_id']}:")
        print(f"  Deforested Area year i: {deforested_area_yeari:.2f} ha/yr")
    
    return total_results
  
#### Step 4A.3: Define equation to calculate Baseline Annual Total Biomass CO2 Emissions ##

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

def calculate_baseline_AD_biomass_co2_emissions(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    
    years = int(scenario_data.get('years', 20))
    if years != 20:
        print(f"Warning: 'years' in the JSON file is set to {years}. Using 20 years for calculations.")
    years = 20
    
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    hist_def_rate = calculate_historical_deforestation_rate(scenario_data)
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = float(subregion["forest_area"]) ### Anthony to check here if forest_area is correct in json and we have a way to get this value
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

#### STEP 4B - Baseline Fire subintervention ############
# if the user selects that intervention_sub_category then the equation for the 
#Baseline fire frequency can be accounted for in the deltaCL in the Baseline, and then that term is 
#not included in the Intervention scenario, assuming there are no more fire disturbances due to the protection.

### Step 4B.1: Define equation to calculate fire disturabance impact using total_carbon_stock obtained in Step 3
def calculate_ldisturbance(total_results, scenario, log_level='info'):
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    results = []
    total_ldisturbance = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
    for subregion in total_results:
        total_carbon_stock = subregion["total_carbon_stock"]
        
        Adisturbance = scenario_data.get("Adisturbance") # area affected by the disturbance
        #Bw = scenario_data.get("Bw") # aboveground biomass of the area affected by the disturance
        #ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        fd = scenario_data.get("fd") # fraction of the biomass affected by the disturbance. using 0.5 for fire in baseline
        
        if None in [Adisturbance, Bw, ratio, fd]:
            if log_level == 'debug':
                print(f"Missing data for calculating Ldisturbance in subregion {aoi_id}. Please check your data.")
            continue
        
        ldisturbance = Adisturbance * total_carbon_stock * fd
        
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

###### Step 4B.2: Define function to multiply baseline fire frequency by ldisturbance
def calculate_fire_co2_emissions(scenario, log_level='info'):
    """
    Calculate fire CO2 emissions based on ldisturbance and fire frequency for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Fire CO2 emissions (tCO2e/yr) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    # Get fire frequency from scenario data
    fire_frequency = float(scenario_data.get("fire_frequency", 0))  # e.g., 0.2 for every 5 years
    
    # Calculate ldisturbance using the existing function
    ldisturbance_results = calculate_ldisturbance(scenario, log_level)
    
    results = []
    total_fire_co2_emissions = 0
    total_area = 0
    
    for ldisturbance_result in ldisturbance_results:
        aoi_id = ldisturbance_result["aoi_id"]
        area = ldisturbance_result["area"]
        ldisturbance = ldisturbance_result["ldisturbance"]
        
        # Calculate fire CO2 emissions
        fire_co2_emissions = ldisturbance * fire_frequency
        
        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "fire_co2_emissions": fire_co2_emissions
        })
        
        total_fire_co2_emissions += fire_co2_emissions
        total_area += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Fire CO2 Emissions: {fire_co2_emissions:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")
    
    average_fire_co2_emissions = total_fire_co2_emissions / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nFire Frequency: {fire_frequency:.4f} (every {1/fire_frequency:.1f} years)")
        print(f"Total Fire CO2 Emissions across all subregions: {total_fire_co2_emissions:.2f} tCO2e/yr")
        print(f"Average Fire CO2 Emissions per hectare: {average_fire_co2_emissions:.2f} tCO2e/ha/yr")
    
    return results

#### STEP 4C - Baseline Fire subintervention ############

## same as equations in Step 5C.1, Step 5C.2 and Step 5C.3 below


#### STEP 5: Define function to estimate average annual biomass growth in Baseline Scenario
def calculate_average_annual_biomass_growth(scenario, log_level='info'):
    """
    Calculate average annual biomass growth (Gtotal) for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Average annual biomass growth (tonnes d.m. ha-1 yr-1) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = float(subregion["area"])
        
        # Retrieve parameters from scenario data
        R = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))   # Root-to-shoot ratio
        #Iv = float(scenario_data.get("Iv", 0))  # Average net annual increment
        #BCEFi = float(scenario_data.get("BCEFi", 0))  # Biomass conversion and expansion factor
        Gw = float(scenario_data.get("Gw", 0)) #average annual increment in t dm/ha/yr
        
        # Calculate Gtotal
        Gtotal = Gw * (1 + R)
        
        subregion_result = {
            "aoi_id": aoi_id,
            "area": area,
            "Gtotal": Gtotal
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area: {area:.2f} ha")
            print(f"  R: {R:.4f}")
            print(f"  Gw: {Gw:.4f} tonnes d.m. ha-1 yr-1")
            print(f"  Gtotal: {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
    
    if log_level == 'debug':
        total_area = sum(result['area'] for result in total_results)
        total_Gtotal = sum(result['area'] * result['Gtotal'] for result in total_results)
        average_Gtotal = total_Gtotal / total_area if total_area > 0 else 0
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Average Gtotal across all subregions: {average_Gtotal:.4f} tonnes d.m. ha-1 yr-1")
    
    return total_results

#### STEP 5.1: Define function to estimate annual increase in biomass carbon stock in Baseline Scenario
def calculate_delta_cg(scenario, log_level='info'):
    """
    Calculate the annual increase in biomass carbon stocks (ΔCG) for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Annual increase in biomass carbon stocks (tC/yr) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    
    # Get the carbon fraction of dry matter from scenario data
    CF = float(scenario_data.get("CF", 0))
    
    # Calculate Gtotal using the existing function
    gtotal_results = calculate_average_annual_biomass_growth(scenario, log_level)
    
    results = []
    total_delta_cg = 0
    total_area = 0
    
    for gtotal_result in gtotal_results:
        aoi_id = gtotal_result["aoi_id"]
        A = gtotal_result["area"]  # Area in hectares
        Gtotal = gtotal_result["Gtotal"]
        
        # Calculate ΔCG
        delta_cg = A * Gtotal * CF * 44/12
        
        results.append({
            "aoi_id": aoi_id,
            "area": A,
            "delta_cg": delta_cg
        })
        
        total_delta_cg += delta_cg
        total_area += A
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area (A): {A:.2f} ha")
            print(f"  Gtotal: {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
            print(f"  Carbon Fraction (CF): {CF:.4f}")
            print(f"  Annual Increase in Biomass Carbon Stocks (ΔCG): {delta_cg:.2f} tCO2e/yr")
    
    average_delta_cg = total_delta_cg / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Total Annual Increase in Biomass Carbon Stocks: {total_delta_cg:.2f} tC/yr")
        print(f"Average Annual Increase in Biomass Carbon Stocks per hectare: {average_delta_cg:.4f} tCO2e/ha/yr")
    
    return results


############### Sub-intervention Calculations ###############
## 5A Protection from Deforestation
## 5B Protection from Fire
## 5C Protection from Illegal Logging

#### STEP 5A: Avoided Deforestation Total Biomass CO2 Calculations #######

###Define equations needed to calculate Avoided Deforestation Emissions from Trees:

## Step 5A.1: Equation to calculate Actual deforested area and Forest area end of year n
def calculate_actual_area_deforested_year_n(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    forest_area = scenario_data["forest_area"] ### Anthony to check here if forest_area is correct in json and we have a way to get this value
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

### Step 5A.2: Calculate area of avoided deforestation in year n
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
    print(f"forest_area_end_of_year: {forest_area_end_of_year:.2f} ha/yr")
    
    return total_results

### Step 5A.3: Calculate Avoided Deforestation Emissions from Trees
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

###Step 5A.4: Define function to estimate average annual biomass growth due to Avoided Deforestation
def calculate_average_annual_biomass_growth(scenario, log_level='info'):
    """
    Calculate average annual biomass growth (Gtotal) for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Average annual biomass growth (tonnes d.m. ha-1 yr-1) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = float(subregion["area"])
        
        # Retrieve parameters from scenario data
        R = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))   # Root-to-shoot ratio
        #Iv = float(scenario_data.get("Iv", 0))  # Average net annual increment
        #BCEFi = float(scenario_data.get("BCEFi", 0))  # Biomass conversion and expansion factor
        Gw = float(scenario_data.get("Gw", 0)) #average annual increment in t dm/ha/yr
        
        # Calculate Gtotal
        Gtotal = Gw * (1 + R)
        
        subregion_result = {
            "aoi_id": aoi_id,
            "area": area,
            "Gtotal": Gtotal
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area: {area:.2f} ha")
            print(f"  R: {R:.4f}")
            print(f"  Gw: {Gw:.4f} tonnes d.m. ha-1 yr-1")
            print(f"  Gtotal: {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
    
    if log_level == 'debug':
        total_area = sum(result['area'] for result in total_results)
        total_Gtotal = sum(result['area'] * result['Gtotal'] for result in total_results)
        average_Gtotal = total_Gtotal / total_area if total_area > 0 else 0
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Average Gtotal across all subregions: {average_Gtotal:.4f} tonnes d.m. ha-1 yr-1")
    
    return total_results

###Step 5A.5: Define function to estimate annual increase in biomass carbon stock due to Avoided Deforestation
def calculate_delta_cg(scenario, log_level='info'):
    """
    Calculate the annual increase in biomass carbon stocks (ΔCG) for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Annual increase in biomass carbon stocks (tC/yr) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    
    # Get the carbon fraction of dry matter from scenario data
    CF = float(scenario_data.get("CF", 0))
    
    # Calculate Gtotal using the existing function
    gtotal_results = calculate_average_annual_biomass_growth(scenario, log_level)
    
    results = []
    total_delta_cg = 0
    total_area = 0
    
    for gtotal_result in gtotal_results:
        aoi_id = gtotal_result["aoi_id"]
        A = gtotal_result["area"]  # Area in hectares
        Gtotal = gtotal_result["Gtotal"]
        
        # Calculate ΔCG
        delta_cg = A * Gtotal * CF
        
        results.append({
            "aoi_id": aoi_id,
            "area": A,
            "delta_cg": delta_cg
        })
        
        total_delta_cg += delta_cg
        total_area += A
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area (A): {A:.2f} ha")
            print(f"  Gtotal: {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
            print(f"  Carbon Fraction (CF): {CF:.4f}")
            print(f"  Annual Increase in Biomass Carbon Stocks (ΔCG): {delta_cg:.2f} tC/yr")
    
    average_delta_cg = total_delta_cg / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Total Annual Increase in Biomass Carbon Stocks: {total_delta_cg:.2f} tC/yr")
        print(f"Average Annual Increase in Biomass Carbon Stocks per hectare: {average_delta_cg:.4f} tC/ha/yr")
    
    return results

#### STEP 5B: Avoided Fire Total Biomass CO2 Calculations #######

## Estimate ∆CL using Equation 2.11: ∆CL = Ldisturbance
# Ldisturbance = annual biomass carbon losses due to disturbances, tonnes C yr-1 (See Equation 2.14)

## Step 5B.1: Define equation to calculate fire disturabance impact using total_carbon_stock obtained in Step 3
def calculate_ldisturbance(total_results, scenario, log_level='info'):
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    results = []
    total_ldisturbance = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
    for subregion in total_results:
        total_carbon_stock = subregion["total_carbon_stock"]
        
        Adisturbance = scenario_data.get("Adisturbance") # area affected by the disturbance
        #Bw = scenario_data.get("Bw") # aboveground biomass of the area affected by the disturance
        #ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        fd = scenario_data.get("fd") # fraction of the biomass affected by the disturbance. using 0.5 for fire in baseline
        
        if None in [Adisturbance, Bw, ratio, fd]:
            if log_level == 'debug':
                print(f"Missing data for calculating Ldisturbance in subregion {aoi_id}. Please check your data.")
            continue
        
        ldisturbance = Adisturbance * total_carbon_stock * fd
        
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

### Step 5B.2: Estimate avoided biomass burning emissions from avoided fire
def calculate_avoided_biomass_burning_emissions(scenario, log_level='info'):
    """
    Calculate avoided biomass burning emissions for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Avoided biomass burning emissions (tCO2/yr) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    Comb_fac = 0.34  # Combustion factor, Value obtained from winrock AFOLU tool, from IPCC (2006)
    bio_ef = float(scenario_data.get("bio_ef", 0))  # Fire emissions factor
    fire_incidence_rate = float(scenario_data.get("fire_incidence_rate", 0))
    
    total_results = []
    total_avoided_emissions = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = float(subregion["area"])  # Forest area in hectares
        tree_carbon_stock = float(subregion.get("tree_carbon_stock", 0))  # Tree Carbon Stock in tC/ha
        
        # Calculate avoided biomass burning emissions
        avoided_emissions = (
            forest_area * 
            (fire_incidence_rate * tree_carbon_stock * (1/0.47) * Comb_fac * bio_ef) * 
            1e-3  # Convert to tCO2/yr
        )
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": forest_area,
            "avoided_biomass_burning_emissions": avoided_emissions
        })
        
        total_avoided_emissions += avoided_emissions
        total_area += forest_area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Forest Area: {forest_area:.2f} ha")
            print(f"  Tree Carbon Stock: {tree_carbon_stock:.2f} tC/ha")
            print(f"  Avoided Biomass Burning Emissions: {avoided_emissions:.4f} tCO2/yr")
    
    average_avoided_emissions = total_avoided_emissions / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nFire Incidence Rate: {fire_incidence_rate:.4f}")
        print(f"Combustion Factor (Comb_fac): {Comb_fac:.2f}")
        print(f"Fire Emissions Factor (bio_ef): {bio_ef:.4f}")
        print(f"Total area across all subregions: {total_area:.2f} ha")
        print(f"Total Avoided Biomass Burning Emissions: {total_avoided_emissions:.4f} tCO2/yr")
        print(f"Average Avoided Biomass Burning Emissions per hectare: {average_avoided_emissions:.4f} tCO2/ha/yr")
    
    return total_results

### Step 5B.3: Define function to calculate biomass growth due to avoided fire
def calculate_average_annual_biomass_growth(scenario, log_level='info'):
    """
    Calculate average annual biomass growth (Gtotal) for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Average annual biomass growth (tonnes d.m. ha-1 yr-1) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = float(subregion["area"])
        
        # Retrieve parameters from scenario data
        R = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))   # Root-to-shoot ratio
        #Iv = float(scenario_data.get("Iv", 0))  # Average net annual increment
        #BCEFi = float(scenario_data.get("BCEFi", 0))  # Biomass conversion and expansion factor
        Gw = float(scenario_data.get("Gw", 0)) #average annual increment in t dm/ha/yr
        
        # Calculate Gtotal
        Gtotal = Gw * (1 + R)
        
        subregion_result = {
            "aoi_id": aoi_id,
            "area": area,
            "Gtotal": Gtotal
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area: {area:.2f} ha")
            print(f"  R: {R:.4f}")
            print(f"  Gw: {Gw:.4f} tonnes d.m. ha-1 yr-1")
            print(f"  Gtotal: {Gtotal:.4f} tonnes d.m. ha-1 yr-1")
    
    if log_level == 'debug':
        total_area = sum(result['area'] for result in total_results)
        total_Gtotal = sum(result['area'] * result['Gtotal'] for result in total_results)
        average_Gtotal = total_Gtotal / total_area if total_area > 0 else 0
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Average Gtotal across all subregions: {average_Gtotal:.4f} tonnes d.m. ha-1 yr-1")
    
    return total_results


#### STEP 5C: Avoided Illegal Logging Total Biomass CO2 Calculations #######

######## Anthony to check how to obtain illegal logging rate (illegal_logging_rate) - it is in the json already

### Step 5C.1: Emissions from Incidental Damage from Illegal Logging
## Define equation to calculate Carbon Emissions from Incidental Damage due to Illegal Logging
def calculate_incidental_damage_illegal(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data["aoi_subregions"]
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        LDF_factor_1 = float(scenario_data.get("LDF_factor_1", 0))
        LDF_factor_2 = float(scenario_data.get("LDF_factor_2", 0))
        illegal_logging_rate = float(scenario_data.get("illegal_logging_rate", 0)) # in json, Anthony to check the source for the illegal_logging_rate
        average_total_tco2e = float(scenario_data.get("average_total_tco2e", 0))
        
        LDF = LDF_factor_1 * average_total_tco2e + LDF_factor_2
        IncidentalDamage_illegal = LDF * illegal_logging_rate
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "LDF": LDF,
            "IncidentalDamage_illegal": IncidentalDamage_illegal
        })
    
    return total_results

### Step 5C.2: Define equation to calculate Illegal Logging Emissions from Trees
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
        TimberTree_illegal = area * illegal_logging_rate * ELE
        
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

### Step 5C.3: Estimate emissions from Community offtake
def calculate_community_offtake_emissions(json_file_path, scenario, log_level='info'):
    data = load_json_data(json_file_path)
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = float(subregion.get("area", 0)) ## Anthony to check 
        community_offtake = float(subregion.get("community_offtake", 0)) ## TO INCLUDE IN THE JSON FILE
        
        D = float(scenario_data.get("D", 0))  # Wood Density
        ElE_factor_1 = float(scenario_data.get("ElE_factor_1", 0))
        ElE_factor_2 = float(scenario_data.get("ElE_factor_2", 0))
        
        # Calculate ELE
        ELE = (ElE_factor_1 * D) - ElE_factor_2
        
        # Calculate TimberTreecommunity
        TimberTreecommunity = forest_area * community_offtake * ELE
        
        # Calculate Incidental DamageComm (assumed to be a fraction of TimberTreecommunity)
        # You may need to adjust this calculation based on your specific requirements
        incidental_damage_fraction = float(scenario_data.get("incidental_damage_fraction", 0)) # added this to the json
        IncidentalDamageComm = TimberTreecommunity * incidental_damage_fraction ## need to check if it's in the json
        
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



####### Step 6: Calculate all - equations for baseline ######
## 6A Baseline Avoided Deforestation
## 6B Baseline Avoided Fire
## 6C Baseline Avoided Illegal Logging

### Step 6A: Calculate all - Baseline Deforestation Scenario #####
def calculate_all_bauAD(scenario, log_level='info'):
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

    annual_biomass_growth_baseline = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["annual_biomass_growth_baseline"] = annual_biomass_growth_baseline
    
    annual_biomass_carbon_stock_gain_baseline = calculate_delta_cg(scenario=scenario, log_level=log_level)
    result["annual_biomass_carbon_stock_gain_baseline"] = annual_biomass_carbon_stock_gain_baseline
    
    annual_baseline_AD_biomass_co2_emissions = calculate_baseline_AD_biomass_co2_emissions(historical_def_rate, carbon_stock_loss, scenario=scenario, log_level=log_level)
    result["annual_baseline_AD_biomass_co2_emissions"] = annual_baseline_AD_biomass_co2_emissions
    if log_level == 'debug':
        print(f"Baseline AD Annual Biomass CO2 emissions: {annual_baseline_AD_biomass_co2_emissions} tCO2e/yr")
    return result

### Step 6B: Calculate all - Baseline Emissions Fire Scenario ####
def calculate_all_bau_AF(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    annual_biomass_growth_baseline = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["annual_biomass_growth_baseline"] = annual_biomass_growth_baseline
    
    annual_biomass_carbon_stock_gain_baseline = calculate_delta_cg(scenario=scenario, log_level=log_level)
    result["annual_biomass_carbon_stock_gain_baseline"] = annual_biomass_carbon_stock_gain_baseline
    
    l_disturbance_baseline = calculate_ldisturbance(total_results, scenario=scenario, log_level=log_level)
    result["l_disturbance_baseline"] = l_disturbance_baseline
    
    avoided_fire_emissions = calculate_avoided_biomass_burning_emissions(scenario=scenario, log_level=log_level)
    result["avoided_fire_emissions"] = avoided_fire_emissions
    if log_level == 'debug':
        print(f"Baseline AF Annual Biomass CO2 emissions: {avoided_fire_emissions} tCO2e/yr")
    return result

### Step 6B: Calculate all - Baseline Emissions Illegal Logging ####

def calculate_all_bau_AIL(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    annual_biomass_growth_baseline = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["annual_biomass_growth_baseline"] = annual_biomass_growth_baseline
    
    annual_biomass_carbon_stock_gain_baseline = calculate_delta_cg(scenario=scenario, log_level=log_level)
    result["annual_biomass_carbon_stock_gain_baseline"] = annual_biomass_carbon_stock_gain_baseline
    
    loss_incidental_damage_illegal = calculate_incidental_damage_illegal(scenario=scenario, log_level=log_level)
    result["incidental_damage_illegal"] = incidental_damage_illegal
    
    loss_illegal_logging_emissions = calculate_logging_emissions_illegal(scenario=scenario, log_level=log_level)
    result["loss_illegal_logging_emissions"] = loss_illegal_logging_emissions
    
    community_offtake_emissions = calculate_community_offtake_emissions(scenario=scenario, log_level=log_level)
    result["community_offtake_emissions"] = community_offtake_emissions
    if log_level == 'debug':
        print(f"Baseline AF Annual Biomass CO2 emissions: {avoided_fire_emissions} tCO2e/yr")
    return result


####### Step 7: Calculate all - equations for interventions ######
## 7A Sub-intervention Avoided Deforestation scenario
## 7B Sub-intervention Avoided Fire scenario
## 7C Sub-intervention Avoided Illegal Logging scenario

#### Step 7A: Calculate all - Avoided Deforestation Scenario 
def calculate_all_int_AD(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    actual_area_deforested_each_year = calculate_actual_area_deforested_year_n(scenario=scenario, log_level=log_level)
    result["actual_area_deforested_each_year"] = actual_area_deforested_each_year
    
    area_avoided_deforestation_year_n = calculate_area_avoided_deforestation_year_n(scenario=scenario, log_level=log_level)
    result["area_avoided_deforestation_year_n"] = area_avoided_deforestation_year_n
    
    yearly_avoided_emissions_trees = calculate_avoided_emissions_trees_each_year(scenario=scenario, log_level=log_level)
    result["yearly_avoided_emissions_trees"] = yearly_avoided_emissions_trees
    
    annual_biomass_growth_AD = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["annual_biomass_growth_AD"] = annual_biomass_growth_AD
    
    annual_biomass_carbon_stock_gain_AD = calculate_delta_cg(scenario=scenario, log_level=log_level)
    result["annual_biomass_carbon_stock_gain_AD"] = annual_biomass_carbon_stock_gain_AD

    yearly_avoided_burning = calculate_avoided_biomass_burning_emissions(scenario=scenario, log_level=log_level)
    result["yearly_decrease_biomass_stock"] = yearly_decrease_biomass_stock
    
    loss_incidental_damage_illegal = calculate_incidental_damage_illegal(scenario=scenario, log_level=log_level)
    result["incidental_damage_illegal"] = incidental_damage_illegal
    
    loss_illegal_logging_emissions = calculate_logging_emissions_illegal(scenario=scenario, log_level=log_level)
    result["loss_illegal_logging_emissions"] = loss_illegal_logging_emissions
    
    community_offtake_emissions = calculate_community_offtake_emissions(scenario=scenario, log_level=log_level)
    result["community_offtake_emissions"] = community_offtake_emissions
    
    return result

#### Step 7B: Calculate all - Avoided Fire Scenario 
def calculate_all_int_AD(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    annual_biomass_growth_AF = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["annual_biomass_growth_AF"] = annual_biomass_growth_AF
    
    annual_biomass_carbon_stock_gain_AF = calculate_delta_cg(scenario=scenario, log_level=log_level)
    result["annual_biomass_carbon_stock_gain_AF"] = annual_biomass_carbon_stock_gain_AF

    yearly_avoided_burning = calculate_avoided_biomass_burning_emissions(scenario=scenario, log_level=log_level)
    result["yearly_decrease_biomass_stock"] = yearly_decrease_biomass_stock
    
    loss_incidental_damage_illegal = calculate_incidental_damage_illegal(scenario=scenario, log_level=log_level)
    result["incidental_damage_illegal"] = incidental_damage_illegal
    
    loss_illegal_logging_emissions = calculate_logging_emissions_illegal(scenario=scenario, log_level=log_level)
    result["loss_illegal_logging_emissions"] = loss_illegal_logging_emissions
    
    community_offtake_emissions = calculate_community_offtake_emissions(scenario=scenario, log_level=log_level)
    result["community_offtake_emissions"] = community_offtake_emissions
    
    return result

#### Step 7C: Calculate all - Avoided Illegal Logging and Community Offtake Scenario 
def calculate_all_int_AD(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    annual_biomass_growth_AIL = calculate_average_annual_biomass_growth(scenario=scenario, log_level=log_level)
    result["annual_biomass_growth_AIL"] = annual_biomass_growth_AIL
    
    annual_biomass_carbon_stock_gain_AIL = calculate_delta_cg(scenario=scenario, log_level=log_level)
    result["annual_biomass_carbon_stock_gain_AIL"] = annual_biomass_carbon_stock_gain_AIL

    loss_incidental_damage_illegal = calculate_incidental_damage_illegal(scenario=scenario, log_level=log_level)
    result["incidental_damage_illegal"] = incidental_damage_illegal
    
    loss_illegal_logging_emissions = calculate_logging_emissions_illegal(scenario=scenario, log_level=log_level)
    result["loss_illegal_logging_emissions"] = loss_illegal_logging_emissions
    
    community_offtake_emissions = calculate_community_offtake_emissions(scenario=scenario, log_level=log_level)
    result["community_offtake_emissions"] = community_offtake_emissions
    
    return result


#### Step 8: Calculate differnce between bau and intervention scenarios #####

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

        # # N2O emissions calculation
        # if "nutrient management" in intervention_subcategory:
        #     FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
        #            temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100)
        #     FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
        #            temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100)
        # else:
        #     FSN = FON = 0
        
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
  
### Step 10.4: Define function to generate output json file

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
