# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: ARR_final.py
# Author: Barbara Bomfim
# Date Started: 06/28/2024
# Last Edited: 17/09/2024
# Purpose: AFOLU GHG Calculations for Afforestation/Reforestation (A/R) Interventions
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data from Dynamic World land use data,  data layers
# User inputs: User inputs depend on sub-intervention category (planted or natural forest; 
#  tillage; nutrient management; fire management). 
# -Baseline data:
#           -land area: area under business-as-usual land use
#           -land_use_type: previous land use type (Agroforestry, Annual cropland, Annual fallow, Degraded land, Flooded rice, Grassland)
#           -forest_type: selection from dropdown list
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#           -fire_used: yes, no
#           -years: intervention durationf (number of years)
#  -All sub-interventions require:
#           -forest_area: area that is being changed from previous land use to new A/R land use
#           -number of improvements over initial scenario (0, 1 or 2)
#.    -For Forest creation sub-intervention:
#           -forest_creation: planted or natural forest (native forest or mangrove forest)
#           -forest_type: selection from dropdown list
#           -H: annual wood removal - merchantable round wood over bark - (m3/yr)
#           -FGtrees = annual volume of fuelwood removal of whole trees (m3/yr)
#           -FGpart = annual volume of fuelwood removal of tree parts (m3/yr)
#     -For tillage sub-intervention:
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#.    -For nutrient management sub-intervention:
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#     -For fire management, user must provide:
#.          -fire_used: yes, no
#           -frequency of burn (number of years)

# Calculations requires the following parameters and their associated uncertainty:

##AGB
#AGBREF:AGB for the climate zone and soil type (tC/ha)

##BGB
#BGBREF:BGB for the climate zone and soil type (tC/ha)

##SOC
# SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
# FLU: Land use emissions factor (1 for planted forest, 1 for natural forest)
# FMG: Management factor for both business as usual and intervention scenario 
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

######### Annual GHG impact calculations ############

## Import libraries
from pygrowthmodels import chapman_richards

import json
import os
from typing import Dict, List, Any # edited this line

import pandas as pd
import numpy as np
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor
import ee
import re
pd.__version__

#Adjusted helpers_ARR file name
from helpers_ARR import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock, sanitize_filename

## Step 0: troubleshooting
# Define function to check for required global variables
def check_global_variables():
    required_globals = ['data', 'polygon', 'average_agbd_dict']
    for var in required_globals:
        if var not in globals():
            raise NameError(f"Required global variable '{var}' is not defined")
          
check_global_variables()

##### STEP 1 - Estimate mean annual AGB (in Mg/ha) using GEDI dataset ######
# GEE path: https://code.earthengine.google.com/?scriptPath=users%2Fbabomfimf%2FAGB-GEDI-L4A%3Atest-3_AGB_annual_mean
# Anthony's GEE: https://code.earthengine.google.com/44d6aaa5db4764b5b5f3825baf900d04

# Authenticate GEE #
# if this doesn't work authenticate from the command line  by running `earthengine authenticate`
ee.Authenticate()

# Initialize Earth Engine
# Use the GEE project name here
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

# Access the GEDI L4A Monthly dataset
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
    scale=25  # Scale should match the resolution of the dataset
)
average_agbd_dict = average_agbd.getInfo()
print(f"{average_agbd_dict} Mg/ha")

print("\nArea of the polygon:")
print(f"{polygon.area().getInfo()/10000} ha") # to get area in hectares

## JSON data input #####

# JSON file path - Check the current working directory
current_directory = os.getcwd()
print("Current working directory:", current_directory)

# Define the relative file path
file_name = 'ARR_final.json'
file_path = os.path.join(current_directory, file_name)
json_file_path = './ARR_final.json'

# Load data from the JSON file
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
# here is the sum of aboveground and belowground carbon stocks
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

#### STEP 4 - Calculate Annual CO2 impact (∆CG) #####

##Annual-CO2-impact = [[Carbon-i - Carbon-0] - [Carbon-bau - DeltaCarbon-0]]/years
#Carbon-0: initial carbon stock (tC/ha) is the average_agbd_cstock calculated above for year 0
#Carbon-i: carbon stock after intervention (tC/ha) - average_agbd_cstock for year i after intervention
#Carbon-bau: carbon stock under business-as-usual (without intervention) - average_agbd_cstock for year i without intervention

## Step 4.1: Richards-Chapman growth model to obtain Carbon-i:
# Define the Chapman-Richards AGB function using pygrowthmodels
# as per https://github.com/drodriguezperez/pygrowthmodels/blob/master/pygrowthmodels/chapman_richards.py#L26

def chapman_richards(time, alpha, beta, rate, slope):
    """
    Parameters
    ----------
    time : time
    alpha : upper asymptote
    beta : growth range
    rate : growth rate
    slope : slope of growth

    Returns:
    float: Aboveground biomass (t d.m./ha).
    """
    return chapman_richards(time, alpha, beta, rate, slope)
  
## Step 4.2: Richards-Chapman growth model to obtain Carbon-i and Carbon-bau (depending on bau scenario):
# Apply Richards-Chapman growth model to calculate AGB annual growth after X number of years

def calculate_agb_growth(scenario, years, log_level='info'):
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    subregions = scenario_data.get("aoi_subregions", [])
 
    # mappings based on:https://github.com/drodriguezperez/pygrowthmodels/blob/master/pygrowthmodels/chapman_richards.py#L26)
    alpha = scenario_data.get("MAX", 0)
    rate = scenario_data.get("k", 0)
    slope = scenario_data.get("m", 0)
    time = scenario_data.get("years", 0) 
    CF = scenario_data.get("CF", 0) 
    beta = 0
    
    results = []
    total_agb_rc = 0
    total_area_rc = 0

    for subregion in subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        agb_rc = chapman_richards(time, alpha, beta, rate, slope)
        
        average_agb_rc_cstock = convert_to_c(agb_rc, CF)
        average_agb_rc_tco2e = convert_to_co2e(average_agb_rc_cstock)
        
        results.append({
            "aoi_id": aoi_id,
            "age": years,
            "area": area,
            "agb_rc": average_agb_rc_tco2e
        })
        
        total_agb_rc += average_agb_rc_tco2e * area
        total_area_rc += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Aboveground Biomass (AGB): {average_agb_rc_tco2e:.2f} tCO2e/ha at age {years}")
    
    average_agb_rc = total_agb_rc / total_area_rc if total_area_rc > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal aboveground biomass across all subregions: {total_agb_rc:.2f} tCO2e")
        print(f"Average aboveground biomass per hectare: {average_agb_rc_tco2e:.2f} tCO2e/ha")
    
    return results

## Step 4.3: Define equation to calculate Annual CO2 impact (∆CG)
# Carbon-i and Carbon-bau are obtained using previous richards-chapman model output

def calculate_annual_co2_impact(scenario, years, average_agbd_tco2e, log_level='info'):
    global data

    # use this one for carbon_i and carbon_bau
    calculate_agb_growth(scenario=scenario, years=years)
    
    carbon_i = get_carbon_stock(data, scenario=scenario, carbon_time_id="carbon_i")
    carbon_0 = average_agbd_tco2e
    carbon_bau = get_carbon_stock(data, scenario=scenario, carbon_time_id="carbon_bau")
    
    results = []
    total_impact = 0
    total_area = 0
    
    for i, subregion in enumerate(carbon_i):
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        carbon_i_co2e = convert_to_co2e(subregion["carbon_stock"]) # from output of previous function
        # retrive based on aoi_id subregion level
        carbon_0_co2e = next((item for item in carbon_0 if item['aoi_id'] == aoi_id), None)["carbon_stock"]
        carbon_bau_co2e = convert_to_co2e(carbon_bau[i]["carbon_stock"]) # from output of previous function
        
        delta_CG = ((carbon_i_co2e - carbon_0_co2e) - (carbon_bau_co2e - carbon_0_co2e)) / years
        
        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "delta_CG": delta_CG
        })
        
        total_impact += delta_CG * area
        total_area += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Annual CO2 impact (∆CG): {delta_CG:.2f} tCO2e/yr in {years} years of intervention")
    
    average_impact = total_impact / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal annual CO2 impact across all subregions: {total_impact:.2f} tCO2e/yr")
        print(f"Average annual CO2 impact per hectare: {average_impact:.2f} tCO2e/ha/yr")
    
    return results

### Step 5 #### Estimate Initial change in biomass carbon stocks on converted land
## Adapting IPCC Equation 2.16 to estimate ∆C-CONVERSION as the immediate change in biomass due to conversion.
## ∆C-CONVERSION = (AGB month after - AGB before) x area of land converted

def estimate_co2_conversion(scenario, log_level='info'):
    """
    Estimate ∆C-CONVERSION as the immediate change in biomass due to conversion for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Immediate change in total biomass carbon due to conversion (tCO2e/yr) for each subregion.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    results = []
    total_delta_co2_conversion = 0
    total_area = 0
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        agb_month_after = scenario_data.get("agb_month_after")
        agb_month_after_co2e = convert_to_co2e(agb_month_after)
        
        agb_before = scenario_data.get("agb_before")
        agb_before_co2e = convert_to_co2e(agb_before)
        
        ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        
        delta_co2_conversion = ((agb_month_after_co2e + (agb_month_after_co2e * ratio)) - 
                                (agb_before_co2e + (agb_before_co2e * ratio))) * area
        
        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "delta_co2_conversion": delta_co2_conversion
        })
        
        total_delta_co2_conversion += delta_co2_conversion
        total_area += area
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Immediate change in biomass carbon: {delta_co2_conversion:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")
    
    average_delta_co2_conversion = total_delta_co2_conversion / total_area if total_area > 0 else 0
    
    if log_level == 'debug':
        print(f"\nTotal immediate change in biomass carbon across all subregions: {total_delta_co2_conversion:.2f} tCO2e/yr")
        print(f"Average immediate change in biomass carbon per hectare: {average_delta_co2_conversion:.2f} tCO2e/ha/yr")
    
    return results    

### Step 6 #### Calculate Annual Decrease in Biomass Stock ####

## Estimate ∆CL using Equation 2.11: ∆CL = Lwood −removals + Lfuelwood + Ldisturbance
# Lwood-removals = annual carbon loss due to wood removals, tonnes C yr-1 (Equation 2.12) 
# Lfuelwood = annual biomass carbon loss due to fuelwood removals, tC/yr (Equation 2.13)
# Ldisturbance = annual biomass carbon losses due to disturbances, tonnes C yr-1 (See Equation 2.14)

## Step 6.1: Define Lwood-removals equation
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
        BCEFr = scenario_data.get("BCEFr")
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

## Step 6.2: Define Lfuelwood equation
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

## Step 6.3: Define Ldisturbance equation
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

## Step 6.4: Define ∆CL equation (∆CL = Lwood −removals + Lfuelwood + Ldisturbance)
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

### Step 7 - Calculate annual change in carbon stocks in biomass on land converted to other land-use category ####

## Estimate ∆CB using IPCC Equation 2.15: ∆CB = ∆CG + ∆C-CONVERSION − ∆CL
#∆CB = annual change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr)
#∆CG = annual increase in carbon stocks in biomass due to growth on land converted to another land-use category (tCO2e/yr)
#∆C-CONVERSION = initial change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr)
#∆CL = annual decrease in biomass carbon stocks due to losses from harvesting, fuel wood gathering and disturbances on land converted to other land-use category (tCO2e/yr)

## Step 7.1: Define ∆CB Equation
#NOTE: if you want to see the prints of these functions use log_level='debug'
# def calculate_annual_change_in_carbon_stocks(delta_CG_results, delta_co2_conversion_results, delta_cl_results, log_level='debug'):
def calculate_annual_change_in_carbon_stocks(delta_CG_results, delta_co2_conversion_results, delta_cl_results, log_level='info'):
    """
    Calculate the annual change in carbon stocks in biomass on land converted to other land-use category for each subregion.
    Parameters:
    delta_CG_results (list): List of dictionaries containing ∆CG results for each subregion.
    delta_co2_conversion_results (list): List of dictionaries containing ∆CO2-conversion results for each subregion.
    delta_cl_results (list): List of dictionaries containing ∆CL results for each subregion.
    log_level (str): Logging level. Set to 'debug' for detailed output.
    Returns:
    list: Annual change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr) for each subregion.
    """
    results = []
    total_annual_change = 0
    total_area = 0

    # Ensure all input lists have the same length
    print(delta_CG_results, delta_co2_conversion_results, delta_cl_results)
    if not (len(delta_CG_results) == len(delta_co2_conversion_results) == len(delta_cl_results)):
        if log_level == 'debug':
            print("Error: Input lists have different lengths. Please ensure all subregions are represented in each input.")
        return None

    for delta_CG, delta_co2_conversion, delta_cl in zip(delta_CG_results, delta_co2_conversion_results, delta_cl_results):
        # Ensure we're dealing with the same subregion in all inputs
        if not (delta_CG['aoi_id'] == delta_co2_conversion['aoi_id'] == delta_cl['aoi_id']):
            if log_level == 'debug':
                print(f"Error: Mismatched subregion IDs: {delta_CG['aoi_id']}, {delta_co2_conversion['aoi_id']}, {delta_cl['aoi_id']}")
            continue

        aoi_id = delta_CG['aoi_id']
        area = delta_CG['area']  # Assuming area is the same in all inputs

        # Calculate annual change for the subregion
        annual_change = delta_CG['delta_CG'] + delta_co2_conversion['delta_co2_conversion'] - delta_cl['delta_cl']

        results.append({
            "aoi_id": aoi_id,
            "area": area,
            "annual_change": annual_change
        })

        total_annual_change += annual_change
        total_area += area

        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Annual change in carbon stocks: {annual_change:.2f} tCO2e/yr")
            print(f"  Area: {area:.2f} ha")

    average_annual_change = total_annual_change / total_area if total_area > 0 else 0

    if log_level == 'debug':
        print(f"\nTotal annual change in carbon stocks across all subregions: {total_annual_change:.2f} tCO2e/yr")
        print(f"Average annual change in carbon stocks per hectare: {average_annual_change:.2f} tCO2e/ha/yr")

    return results

### Step 8: Calculate all ARR variables ####

# ## Define equation to calculate all calculation steps using equations defined above
# def calculate_all(scenario, log_level='info'):
#     result = {}
#     
#     average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
#     result["average_agbd_tco2e"] = average_agbd_tco2e
#     
#     mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
#     result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
#     
#     delta_CG = calculate_annual_co2_impact(scenario=scenario, years=20, average_agbd_tco2e=average_agbd_tco2e, log_level=log_level)
#     result["delta_CG"] = delta_CG
#     
#     delta_co2_conversion = estimate_co2_conversion(scenario=scenario, log_level=log_level)
#     result["delta_co2_conversion"] = delta_co2_conversion
#     
#     # Calculate Lwood-removals
#     lwood_removals = calculate_lwood_removals(scenario=scenario, log_level=log_level)
#     result["lwood_removals"] = lwood_removals
#     if log_level == 'debug':
#         print(f"Annual biomass carbon loss due to wood removals (Lwood-removals): {lwood_removals} tCO2e/yr")
#     
#     # Calculate Lfuelwood
#     lfuelwood = calculate_lfuelwood(scenario=scenario, log_level=log_level)
#     result["lfuelwood"] = lfuelwood
#     if log_level == 'debug':
#         print(f"Annual biomass carbon loss due to fuelwood removals (Lfuelwood): {lfuelwood} tCO2e/yr")
#     
#     # Calculate Ldisturbance
#     ldisturbance = calculate_ldisturbance(scenario=scenario, log_level=log_level)
#     result["ldisturbance"] = ldisturbance
#     if log_level == 'debug':
#         print(f"Annual biomass carbon losses due to disturbances (Ldisturbance): {ldisturbance} tCO2e/yr")
#     
#     # Calculate ∆CL
#     delta_cl = calculate_delta_cl(lwood_removals, lfuelwood, ldisturbance, log_level=log_level)
#     result["delta_cl"] = delta_cl
#     if log_level == 'debug':
#         print(f"Biomass stock Loss (∆CL): {delta_cl} tCO2e/yr")
#     
#     # Calculate ∆CB - the annual change in carbon stocks
#     annual_change_in_carbon_stocks = calculate_annual_change_in_carbon_stocks(delta_CG, delta_co2_conversion, delta_cl, log_level=log_level)
#     result["annual_change_in_carbon_stocks"] = annual_change_in_carbon_stocks
#     if log_level == 'debug':
#         print(f"Annual change in total biomass carbon stock (∆CB): {annual_change_in_carbon_stocks} tCO2e/yr")
#     
#     return result

# biomass = {}
# biomass["business_as_usual"] = calculate_all(scenario="business_as_usual")
# biomass["intervention"] = calculate_all(scenario="intervention")
# 
# biomass_co2_result = {}
# 
# for inter_item, bau_item in zip(biomass["intervention"]["annual_change_in_carbon_stocks"], 
#                                 biomass["business_as_usual"]["annual_change_in_carbon_stocks"]):
#     if inter_item['aoi_id'] != bau_item['aoi_id']:
#         raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
#     
#     aoi_id = inter_item['aoi_id']
#     difference = inter_item['annual_change'] - bau_item['annual_change']
#     
#     biomass_co2_result[aoi_id] = difference# inter_result = calculate_all(scenario="intervention")
# 
# print(biomass_co2_result)
# 
# biomass_co2_result_error_positive = 10  # this is a percentage
# 
# biomass_co2_sd = {}
# 
# for aoi_id, difference in biomass_co2_result.items():
#     # Calculate SD for each aoi_id
#     sd = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96
#     biomass_co2_sd[aoi_id] = sd

# #### edited code to fix issues
# def calculate_all(scenario, log_level='info'):
#     result = {}
#     
#     try:
#         check_global_variables()
#         
#         average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
#         result["average_agbd_tco2e"] = average_agbd_tco2e
#         
#         mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
#         result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
#         
#         delta_CG = calculate_annual_co2_impact(scenario=scenario, years=20, average_agbd_tco2e=average_agbd_tco2e, log_level=log_level)
#         result["delta_CG"] = delta_CG
#         
#         delta_co2_conversion = estimate_co2_conversion(scenario=scenario, log_level=log_level)
#         result["delta_co2_conversion"] = delta_co2_conversion
#         
#         lwood_removals = calculate_lwood_removals(scenario=scenario, log_level=log_level)
#         result["lwood_removals"] = lwood_removals
#         
#         lfuelwood = calculate_lfuelwood(scenario=scenario, log_level=log_level)
#         result["lfuelwood"] = lfuelwood
#         
#         ldisturbance = calculate_ldisturbance(scenario=scenario, log_level=log_level)
#         result["ldisturbance"] = ldisturbance
#         
#         delta_cl = calculate_delta_cl(lwood_removals, lfuelwood, ldisturbance, log_level=log_level)
#         result["delta_cl"] = delta_cl
#         
#         annual_change_in_carbon_stocks = calculate_annual_change_in_carbon_stocks(delta_CG, delta_co2_conversion, delta_cl, log_level=log_level)
#         result["annual_change_in_carbon_stocks"] = annual_change_in_carbon_stocks
#         
#     except Exception as e:
#         print(f"Error in calculate_all for scenario {scenario}: {str(e)}")
#         result["error"] = str(e)
#     
#     return result

#### edited calculate_all function to fix issues
def calculate_all(scenario, log_level='info'):
    result = {}
    
    try:
        check_global_variables()
        
        steps = [
            ('average_agbd_tco2e', lambda: mean_annual_agc_stock(scenario=scenario, log_level=log_level)),
            ('mean_annual_tot_c_stock', lambda: mean_annual_tot_c_stock(result['average_agbd_tco2e'], scenario=scenario, log_level=log_level)),
            ('delta_CG', lambda: calculate_annual_co2_impact(scenario=scenario, years=20, average_agbd_tco2e=result['average_agbd_tco2e'], log_level=log_level)),
            ('delta_co2_conversion', lambda: estimate_co2_conversion(scenario=scenario, log_level=log_level)),
            ('lwood_removals', lambda: calculate_lwood_removals(scenario=scenario, log_level=log_level)),
            ('lfuelwood', lambda: calculate_lfuelwood(scenario=scenario, log_level=log_level)),
            ('ldisturbance', lambda: calculate_ldisturbance(scenario=scenario, log_level=log_level)),
            ('delta_cl', lambda: calculate_delta_cl(result['lwood_removals'], result['lfuelwood'], result['ldisturbance'], log_level=log_level)),
            ('annual_change_in_carbon_stocks', lambda: calculate_annual_change_in_carbon_stocks(result['delta_CG'], result['delta_co2_conversion'], result['delta_cl'], log_level=log_level))
        ]
        
        for key, func in steps:
            try:
                result[key] = func()
                if log_level == 'debug':
                    print(f"Calculated {key} successfully")
            except Exception as step_error:
                print(f"Error calculating {key} for scenario {scenario}: {str(step_error)}")
                result[f"{key}_error"] = str(step_error)
        
    except Exception as e:
        print(f"Unexpected error in calculate_all for scenario {scenario}: {str(e)}")
        result["unexpected_error"] = str(e)
    
    return result

#### Step 9: Calculate bau and intervention scenarios and difference between both

biomass = {
    "business_as_usual": calculate_all(scenario="business_as_usual"),
    "intervention": calculate_all(scenario="intervention")
}

biomass_co2_result = {}
biomass_co2_sd = {}
biomass_co2_result_error_positive = 10  # this is a percentage

# for inter_item, bau_item in zip(biomass["intervention"]["annual_change_in_carbon_stocks"], 
#                                 biomass["business_as_usual"]["annual_change_in_carbon_stocks"]):
#     if inter_item['aoi_id'] != bau_item['aoi_id']:
#         raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
#     
#     aoi_id = inter_item['aoi_id']
#     difference = inter_item['annual_change'] - bau_item['annual_change']
#     
#     biomass_co2_result[aoi_id] = difference
#     biomass_co2_sd[aoi_id] = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96
# 
# print("Biomass CO2 Result:", biomass_co2_result)
# print("Biomass CO2 Standard Deviation:", biomass_co2_sd)

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


#### Step 10: SOIL GHG CALCULATIONS ######

##Step 10.1 Getting JSON data ready to be calculated ##
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

## Step 10.2: Define function to organize results

# def average_results_sims(input_dict):
#     results_unc = {}
#     for key, value in input_dict.items():
#         if isinstance(value, list):
#             if value:  # Check if the list is not empty
#                 average = sum(value) / len(value)
#                 results_unc[key] = average
#             else:
#                 results_unc[key] = 0  # or None, for empty lists
#         else:
#             results_unc[key] = value  # Keep non-list values as they are
#     return results_unc

#edited
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

# def GHGcalc(aoi_id, df, nx, intervention_subcategory, biomass_co2_result):
#     
#     temp_bau = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "business-as-usual")]
#     temp_int = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "intervention")]
#     # Create a structured array that can hold both strings and floats
#     results_unc = {}
#     
#     results_unc['AOI'] = aoi_id
#     results_unc['Area'] = temp_bau['area'].values[0]
#     results_unc['SOC'] = []
#     results_unc['totalC'] = []
#     results_unc[f'N2O_{aoi_id}'] = []
#     # this needs to be calculated with fire management below
#     #CH4 calculation needed
#     # results_unc[f'CH4_{aoi_id}'] = []
#     
#     for m in range(nx):
#         # SOC Stock Change
#         SOCREF = np.random.normal(temp_bau['SOC_ref_tonnes_C_ha'].values[0], temp_bau['SOC_ref_tonnes_C_ha_sd'].values[0])
#         FLUbau = temp_bau['FLU'].values[0]
#         FMGbau = np.random.normal(temp_bau['FMG'].values[0], temp_bau['FMG_sd'].values[0])
#         FIbau = np.random.normal(temp_bau['FI'].values[0], temp_bau['FI_sd'].values[0])
#         FLUint = temp_int['FLU'].values[0]
#         FMGint = np.random.normal(temp_int['FMG'].values[0], temp_int['FMG_sd'].values[0])
#         FIint = np.random.normal(temp_int['FI'].values[0], temp_int['FI_sd'].values[0])
#         SOCbau = SOCREF * FLUbau * FMGbau * FIbau
#         SOCint = SOCREF * FLUint * FMGint * FIint
#         dSOC = SOCbau - SOCint
#         results_unc["SOC"].append(dSOC)
# 
#         # here we add the biomass specific to the aoi_id
#         # NOTE: this number seems too large - check please
#         results_unc["totalC"].append(dSOC * 44/12 + biomass_co2_result[aoi_id]) 
# 
#         # N2O emissions 
#         # Calculate change in N sources
#         FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
#                temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100) if "nutrient management" in intervention_subcategory else 0
#         FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
#                temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100) if "nutrient management" in intervention_subcategory else 0
#         FSOM = dSOC / 10 * 1000  # FSOM
# 
#         # Calculate N emissions
#         EFdir = np.random.normal(temp_bau['ef_direct_n2o'].values[0], temp_bau['ef_direct_n2o_sd'].values[0])
#         FRAC_GASF = np.random.normal(temp_bau['frac_syn_fert'].values[0], temp_bau['frac_syn_fert_sd'].values[0])
#         FRAC_GASM = np.random.normal(temp_bau['frac_org_fert'].values[0], temp_bau['frac_org_fert_sd'].values[0])
#         EF_vol = np.random.normal(temp_bau['ef_vol'].values[0], temp_bau['ef_vol_sd'].values[0])
#         FRAC_LEACH = np.random.normal(temp_bau['frac_leach'].values[0], temp_bau['frac_leach_sd'].values[0])
#         EF_leach = np.random.normal(temp_bau['ef_leaching'].values[0], temp_bau['ef_leaching_sd'].values[0])
#         dirN = (FSN + FON + FSOM) * EFdir  # Direct emissions
#         volN = (FSN * FRAC_GASF + (FON ) * FRAC_GASM) * EF_vol  # Indirect volatilization
#         leachN = (FSN + FON + FSOM) * FRAC_LEACH * EF_leach  # indirect leaching
#         N2O = (dirN + volN + leachN) / 1000 * 44 / 28  # sum and convert to tN2O
# 
#         results_unc[f"N2O_{aoi_id}"].append(N2O)
#         
#         #Edit to include: CH4 calculation
# 
#         # Add change in N2O due to fire management
#         # if "change in fire management" in intervention_subcategory:
#         #     fire_n2o_ef = np.random.normal(temp_bau['burning_n2o_ef_mean'].values[0], temp_bau['burning_n2o_ef_sd'].values[0])
#         #     # Add fires for bau scenario
#         #     if temp_bau['fire_used'].values[0] == "True":
#         #         CF_bau = np.random.normal(temp_bau['combustion_factor_mean'].values[0], temp_bau['combustion_factor_sd'].values[0])
#         #         MB_bau = np.random.normal(temp_bau['fuel_biomass_mean'].values[0], temp_bau['fuel_biomass_sd'].values[0])
#         #         fireN2O_bau = MB_bau * CF_bau * fire_n2o_ef / 1000
#         #         N2O = N2O - fireN2O_bau
#         #         fire_per_bau = temp_bau['fire_management_years'].values[0]
#         #         fire_yrs_bau = [y for y in range(1, 20) if y % fire_per_bau == 0]
#         #         results_unc[m, 4 + np.array(fire_yrs_bau)] = N2O - fireN2O_bau
#         #     if temp_int['fire_used'].values[0] == "True":
#         #         CF_int = np.random.normal(temp_int['combustion_factor_mean'].values[0], temp_int['combustion_factor_sd'].values[0])
#         #         MB_int = np.random.normal(temp_int['fuel_biomass_mean'].values[0], temp_int['fuel_biomass_sd'].values[0])
#         #         fireN2O_int = MB_int * CF_int * fire_n2o_ef / 1000
#         #         results_unc[m, 4] += fireN2O_int
#         #         fire_per_int = temp_int['fire_management_years'].values[0]
#         #         fire_yrs_int = [y for y in range(1, 20) if y % fire_per_int == 0]
#         #         results_unc[m, 4 + np.array(fire_yrs_int)] += fireN2O_int
#         
#         #Edit if incorporating fire management 
#         # # emissions from fire management
#         # if "change in fire management" in intervention_subcategory:
#         #     fire_ch4_ef = np.random.normal(temp_bau['burning_ch4_ef_mean'].values[0], temp_bau['burning_ch4_ef_sd'].values[0])
#         #     # Add fires for bau scenario
#         #     if temp_bau['fire_used'].values[0] == "True":
#         #         fireCH4_bau = MB_bau * CF_bau * fire_ch4_ef / 1000
#         #         results_unc[m, 24] = CH4_live - fireCH4_bau
#         #         results_unc[m, 24 + np.array(fire_yrs_bau)] = CH4_live - fireCH4_bau
#         #     if temp_int['fire_used'].values[0] == "True":
#         #         fireCH4_int = MB_int * CF_int * fire_ch4_ef / 1000
#         #         results_unc[m, 24] += fireCH4_int
#         #         results_unc[m, 24 + np.array(fire_yrs_int)] += fireCH4_int
# 
# 
#     # column_names = ['AOI', 'Area', 'Rep', 'SOC', 'totalC', f'N2O_y{aoi_id}', f'CH4_y{aoi_id}']
#     print("UNC", results_unc)
#     results_unc_df = pd.DataFrame(results_unc)
#     return results_unc_df

####edited code
# def GHGcalc(aoi_id: str, df: pd.DataFrame, nx: int, intervention_subcategory: str, biomass_co2_result: Dict[str, float]) -> pd.DataFrame:
#     temp_bau = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "business-as-usual")]
#     temp_int = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "intervention")]
#     
#     results_unc = {
#         'AOI': aoi_id,
#         'Area': temp_bau['area'].values[0],
#         'SOC': [],
#         'totalC': [],
#         f'N2O_{aoi_id}': [],
#         f'CH4_{aoi_id}': []
#     }
#     
#     for _ in range(nx):
#         # SOC Stock Change calculation
#         SOCREF = np.random.normal(temp_bau['SOC_ref_tonnes_C_ha'].values[0], temp_bau['SOC_ref_tonnes_C_ha_sd'].values[0])
#         FLUbau, FMGbau, FIbau = temp_bau[['FLU', 'FMG', 'FI']].values[0]
#         FLUint, FMGint, FIint = temp_int[['FLU', 'FMG', 'FI']].values[0]
#         
#         SOCbau = SOCREF * FLUbau * np.random.normal(FMGbau, temp_bau['FMG_sd'].values[0]) * np.random.normal(FIbau, temp_bau['FI_sd'].values[0])
#         SOCint = SOCREF * FLUint * np.random.normal(FMGint, temp_int['FMG_sd'].values[0]) * np.random.normal(FIint, temp_int['FI_sd'].values[0])
#         dSOC = SOCbau - SOCint
#         results_unc["SOC"].append(dSOC)
#         
#         # Total carbon calculation
#         results_unc["totalC"].append(dSOC * 44/12 + biomass_co2_result[aoi_id])
# 
#         # N2O emissions calculation
#         if "nutrient management" in intervention_subcategory:
#             FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
#                    temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100)
#             FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
#                    temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100)
#         else:
#             FSN = FON = 0
#         
#         FSOM = dSOC / 10 * 1000
#         
#         EFdir = np.random.normal(temp_bau['ef_direct_n2o'].values[0], temp_bau['ef_direct_n2o_sd'].values[0])
#         FRAC_GASF = np.random.normal(temp_bau['frac_syn_fert'].values[0], temp_bau['frac_syn_fert_sd'].values[0])
#         FRAC_GASM = np.random.normal(temp_bau['frac_org_fert'].values[0], temp_bau['frac_org_fert_sd'].values[0])
#         EF_vol = np.random.normal(temp_bau['ef_vol'].values[0], temp_bau['ef_vol_sd'].values[0])
#         FRAC_LEACH = np.random.normal(temp_bau['frac_leach'].values[0], temp_bau['frac_leach_sd'].values[0])
#         EF_leach = np.random.normal(temp_bau['ef_leaching'].values[0], temp_bau['ef_leaching_sd'].values[0])
#         
#         dirN = (FSN + FON + FSOM) * EFdir
#         volN = (FSN * FRAC_GASF + FON * FRAC_GASM) * EF_vol
#         leachN = (FSN + FON + FSOM) * FRAC_LEACH * EF_leach
#         N2O = (dirN + volN + leachN) / 1000 * 44 / 28
#         
#         results_unc[f"N2O_{aoi_id}"].append(N2O)
# 
#         # CH4 emissions calculation
#         CH4_live = 0  # Initialize CH4 emissions
# 
#         # Fire management calculations
#         if "change in fire management" in intervention_subcategory:
#             fire_n2o_ef = np.random.normal(temp_bau['burning_n2o_ef_mean'].values[0], temp_bau['burning_n2o_ef_sd'].values[0])
#             fire_ch4_ef = np.random.normal(temp_bau['burning_ch4_ef_mean'].values[0], temp_bau['burning_ch4_ef_sd'].values[0])
#             
#             # BAU scenario
#             if temp_bau['fire_used'].values[0] == "True":
#                 CF_bau = np.random.normal(temp_bau['combustion_factor_mean'].values[0], temp_bau['combustion_factor_sd'].values[0])
#                 MB_bau = np.random.normal(temp_bau['fuel_biomass_mean'].values[0], temp_bau['fuel_biomass_sd'].values[0])
#                 fireN2O_bau = MB_bau * CF_bau * fire_n2o_ef / 1000
#                 fireCH4_bau = MB_bau * CF_bau * fire_ch4_ef / 1000
#                 N2O -= fireN2O_bau
#                 CH4_live -= fireCH4_bau
#             
#             # Intervention scenario
#             if temp_int['fire_used'].values[0] == "True":
#                 CF_int = np.random.normal(temp_int['combustion_factor_mean'].values[0], temp_int['combustion_factor_sd'].values[0])
#                 MB_int = np.random.normal(temp_int['fuel_biomass_mean'].values[0], temp_int['fuel_biomass_sd'].values[0])
#                 fireN2O_int = MB_int * CF_int * fire_n2o_ef / 1000
#                 fireCH4_int = MB_int * CF_int * fire_ch4_ef / 1000
#                 N2O += fireN2O_int
#                 CH4_live += fireCH4_int
# 
#         results_unc[f"CH4_{aoi_id}"].append(CH4_live)
# 
#     return pd.DataFrame(results_unc)

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

# def generate_output(input_json):
#     global biomass_co2_result
#     json_data, df, AOIs, intervention_subcategory = load_json_data(input_json)
#
#     # Run Monte Carlo simulations
#     np.random.seed(1)
#     mc = [GHGcalc(aoi_id, df, 1000, intervention_subcategory, biomass_co2_result) for aoi_id in AOIs]
#
#     # Process results
#     result_list = []
#     for iteration in mc:
#         aoi_id = iteration['AOI'].iloc[0]  # Assuming 'AOI' is constant for each iteration
#         result_dict = {'AOI': aoi_id}
#
#         for column in iteration.columns:
#             if column in ['AOI', 'Area', 'Rep']:
#                 result_dict[column] = iteration[column].iloc[0]  # Assuming these are constant for each iteration
#             else:
#                 result_dict[column] = iteration[column].mean()
#                 result_dict[column+'_sd'] = iteration[column].std()
#
#         result_list.append(result_dict)
#
#     # Convert the list of dictionaries to a DataFrame
#     result_df = pd.DataFrame(result_list)
#
#     print(result_df.columns)
#
#     # result_df has the following columns with aoi_id baked in
#     # Run code below to roll them up for totals as well instead of aoi_id
#     # Index(['AOI', 'Area', 'Rep', 'SOC', 'SOC_sd', 'totalC', 'totalC_sd',
#     #    'N2O_tropicalmoist_acrisols', 'N2O_tropicalmoist_acrisols_sd',
#     #    'CH4_tropicalmoist_acrisols', 'CH4_tropicalmoist_acrisols_sd',
#     #    'N2O_tropicalmoist_ferralsols', 'N2O_tropicalmoist_ferralsols_sd',
#     #    'CH4_tropicalmoist_ferralsols', 'CH4_tropicalmoist_ferralsols_sd'],
#     #   dtype='object')
#
#     # Summing specific columns
#     result_df['N2O_thayr'] = result_df['N2O_tropicalmoist_acrisols'] + result_df['N2O_tropicalmoist_ferralsols']
#     result_df['sdN2Oha'] = result_df['N2O_tropicalmoist_acrisols_sd'] + result_df['N2O_tropicalmoist_ferralsols_sd']
#     # result_df['CH4_thayr'] = result_df['CH4_tropicalmoist_acrisols'] + result_df['CH4_tropicalmoist_ferralsols']
#     # result_df['sdCH4ha'] = result_df['CH4_tropicalmoist_acrisols_sd'] + result_df['CH4_tropicalmoist_ferralsols_sd']
#     result_df['CO2_thayr'] = result_df['SOC'] + result_df['totalC']
#     result_df['sdCO2ha'] = result_df['SOC_sd'] + result_df['totalC_sd']
#     result_df['CO2_tyr'] = result_df['CO2_thayr'] * result_df['Area']
#     result_df['sdCO2'] = result_df['sdCO2ha'] * result_df['Area']
#     result_df['N2O_tyr'] = result_df['N2O_thayr'] * result_df['Area']
#     result_df['sdN2O'] = result_df['sdN2Oha'] * result_df['Area']
#     # result_df['CH4_tyr'] = result_df['CH4_thayr'] * result_df['Area']
#     # result_df['sdCH4'] = result_df['sdCH4ha'] * result_df['Area']
#
#     # Drop the individual columns after summing
#     result_df = result_df.drop(columns=[
#         'SOC', 'totalC', 'SOC_sd', 'totalC_sd',
#         'N2O_tropicalmoist_acrisols', 'N2O_tropicalmoist_ferralsols',
#         'N2O_tropicalmoist_acrisols_sd', 'N2O_tropicalmoist_ferralsols_sd',
#         # 'CH4_tropicalmoist_acrisols', 'CH4_tropicalmoist_ferralsols',
#         # 'CH4_tropicalmoist_acrisols_sd', 'CH4_tropicalmoist_ferralsols_sd'
#     ])
#
#     # Aggregate the DataFrame by summing up the columns
#     totals = result_df.sum(numeric_only=True).to_frame().transpose()
#     totals['AOI'] = 'Total'
#     totals['Rep'] = 'Total'
#     result_df = result_df.append(totals, ignore_index=True)
#
#     # Reorder columns to place 'AOI' first
#     cols = ['AOI', 'Rep'] + [col for col in result_df.columns if col not in ['AOI', 'Rep']]
#     result_df = result_df[cols]
#
#     print(result_df)
#     #
#     #NOTE: this is creating all the columns only with projections from year 1 up to you year 20
#     column_names = ['AOI', 'Area', 'CO2_thayr'] + \
#                    [f'N2O_thayr_y{i}' for i in range(1, 21)] + \
#                    ['sdCO2ha'] + [f'sdN2Oha_y{i}' for i in range(1, 21)] + \
#                    ['CO2_tyr'] + [f'N2O_tyr_y{i}' for i in range(1, 21)] + \
#                    ['sdCO2'] + [f'sdN2O_y{i}' for i in range(1, 21)]
#                 # add: CH4 calculation to be added
#                    # [f'CH4_thayr_y{i}' for i in range(1, 21)] + \
#                    # [f'sdCH4ha_y{i}' for i in range(1, 21)] + \
#                    # [f'CH4_tyr_y{i}' for i in range(1, 21)] + \
#                    # [f'sdCH4_y{i}' for i in range(1, 21)]
#
#     resdf = pd.DataFrame(result_df, columns=column_names)
#
#     #NOTE: setting the numbers for projections
#     # # Convert to JSON output
#     restib = []
#     for r in range(len(resdf)):
#         aoi_dict = {
#             'AOI': f'Sub_AOI-{r+1}',
#             'CO2_thayr': [resdf.iloc[r]['CO2_thayr']] * 20,
#             # NOTE: biomass_co2_sd[aoi_id] gives the sd for biomass_co2_result for that aoi_id
#             'sdCO2ha': [resdf.iloc[r]['sdCO2ha']] * 20,
#             'CO2_tyr': [resdf.iloc[r]['CO2_tyr']] * 20,
#             'sdCO2': [resdf.iloc[r]['sdCO2']] * 20,
#             'N2O_thayr': resdf.iloc[r][[f'N2O_thayr_y{i}' for i in range(1, 21)]].tolist(),
#             'sdN2Oha': resdf.iloc[r][[f'sdN2Oha_y{i}' for i in range(1, 21)]].tolist(),
#             'N2O_tyr': resdf.iloc[r][[f'N2O_tyr_y{i}' for i in range(1, 21)]].tolist(),
#             'sdN2O': resdf.iloc[r][[f'sdN2O_y{i}' for i in range(1, 21)]].tolist(),
#                 #TODO: CH4 calculation needed
#             # 'CH4_thayr': resdf.iloc[r][[f'CH4_thayr_y{i}' for i in range(1, 21)]].tolist(),
#             # 'sdCH4ha': resdf.iloc[r][[f'sdCH4ha_y{i}' for i in range(1, 21)]].tolist(),
#             # 'CH4_tyr': resdf.iloc[r][[f'CH4_tyr_y{i}' for i in range(1, 21)]].tolist(),
#             # 'sdCH4': resdf.iloc[r][[f'sdCH4_y{i}' for i in range(1, 21)]].tolist()
#         }
#         restib.append(aoi_dict)
#
#     # Write the DataFrame to a CSV file
#     result_df.to_csv('result.csv', index=False)
#     #
#     json_output = json.dumps(result_df.to_dict(orient='records'), indent=2)
#     #
#     # # # Ensure the output directory exists
#     output_dir = 'output'
#     os.makedirs(output_dir, exist_ok=True)
#     #
#     # # Write the JSON output
#     sanitized_filename = 'output_file.json'  # Replace with actual filename sanitization if needed
#     output_file = os.path.join(output_dir, sanitized_filename)
#     with open(output_file, 'w') as f:
#         json.dump(json_data, f, indent=4)
#
#     sanitized = sanitize_filename(f'{intervention_subcategory}.json')
#     #
#     # # # Write the JSON output
#     output_file = os.path.join(output_dir, sanitized)
#     with open(output_file, 'w') as f:
#        f.write(json_output)
#
#     print(f"\n\nGenerated file: {sanitized_filename} in the output folder.")
#
# generate_output(json_file_path)

# ## improved function to generate output
# def generate_output(input_json: str):
#     global biomass_co2_result
#     json_data, df, AOIs, intervention_subcategory = load_json_data(input_json)
# 
#     np.random.seed(1)
#     mc = [GHGcalc(aoi_id, df, 1000, intervention_subcategory, biomass_co2_result) for aoi_id in AOIs]
#     mc = [result for result in mc if not result.empty] # to filter out empty dataframes
# 
#     if not mc:
#         print{"Error: Not valid results generated. Check input data."}
#         return None
# 
#     result_list = []
#     for iteration in mc:
#         aoi_id = iteration['AOI'].iloc[0]
#         result_dict = {'AOI': aoi_id, 'Area': iteration['Area'].iloc[0]}
# 
#         for column in iteration.columns:
#             if column not in ['AOI', 'Area']:
#                 result_dict[column] = iteration[column].mean()
#                 result_dict[f"{column}_sd"] = iteration[column].std()
# 
#         result_list.append(result_dict)
# 
#     result_df = pd.DataFrame(result_list)
# 
#     # Calculate totals
#     result_df['N2O_thayr'] = result_df[[col for col in result_df.columns if col.startswith('N2O_') and not col.endswith('_sd')]].sum(axis=1)
#     result_df['sdN2Oha'] = np.sqrt((result_df[[col for col in result_df.columns if col.startswith('N2O_') and col.endswith('_sd')]]**2).sum(axis=1))
#     result_df['CH4_thayr'] = result_df[[col for col in result_df.columns if col.startswith('CH4_') and not col.endswith('_sd')]].sum(axis=1)
#     result_df['sdCH4ha'] = np.sqrt((result_df[[col for col in result_df.columns if col.startswith('CH4_') and col.endswith('_sd')]]**2).sum(axis=1))
#     result_df['CO2_thayr'] = result_df['SOC'] + result_df['totalC']
#     result_df['sdCO2ha'] = np.sqrt(result_df['SOC_sd']**2 + result_df['totalC_sd']**2)
#     result_df['CO2_tyr'] = result_df['CO2_thayr'] * result_df['Area']
#     result_df['sdCO2'] = result_df['sdCO2ha'] * result_df['Area']
#     result_df['N2O_tyr'] = result_df['N2O_thayr'] * result_df['Area']
#     result_df['sdN2O'] = result_df['sdN2Oha'] * result_df['Area']
#     result_df['CH4_tyr'] = result_df['CH4_thayr'] * result_df['Area']
#     result_df['sdCH4'] = result_df['sdCH4ha'] * result_df['Area']
# 
#     # Drop individual columns after totaling
#     columns_to_drop = [col for col in result_df.columns if (col.startswith('N2O_') or col.startswith('CH4_')) and col not in ['N2O_thayr', 'CH4_thayr']] + ['SOC', 'totalC', 'SOC_sd', 'totalC_sd']
#     result_df = result_df.drop(columns=columns_to_drop)
# 
#     # Add total row
#     totals = result_df.sum(numeric_only=True).to_frame().T
#     totals['AOI'] = 'Total'
#     result_df = pd.concat([result_df, totals], ignore_index=True)
# 
#     # Reorder columns
#     cols = ['AOI', 'Area'] + [col for col in result_df.columns if col not in ['AOI', 'Area']]
#     result_df = result_df[cols]
# 
#     # Prepare JSON output
#     json_output = []
#     for _, row in result_df.iterrows():
#          aoi_dict = {
#             'AOI': row['AOI'],
#             'Area': row['Area'],
#             'CO2_thayr': [row['CO2_thayr']] * 20,
#             'sdCO2ha': [row['sdCO2ha']] * 20,
#             'CO2_tyr': [row['CO2_tyr']] * 20,
#             'sdCO2': [row['sdCO2']] * 20,
#             'N2O_thayr': [row['N2O_thayr']] * 20,
#             'sdN2Oha': [row['sdN2Oha']] * 20,
#             'N2O_tyr': [row['N2O_tyr']] * 20,
#             'sdN2O': [row['sdN2O']] * 20,
#             'CH4_thayr': [row['CH4_thayr']] * 20,
#             'sdCH4ha': [row['sdCH4ha']] * 20,
#             'CH4_tyr': [row['CH4_tyr']] * 20,
#             'sdCH4': [row['sdCH4']] * 20
#         }
#         json_output.append(aoi_dict)
# 
#     # Generate output files
#      output_dir = 'output'
#      os.makedirs(output_dir, exist_ok=True)
# 
#     # Write CSV output (optional)
#      csv_file = os.path.join(output_dir, 'result.csv')
#      #result_df.to_csv(os.path.join(output_dir, 'result.csv'), index=False)
#      result_df.to_csv(csv_file, index=False)
#      print(f"Generated CSV file: {csv_file}")
# 
#     # Write JSON output
#      sanitized_filename = sanitize_filename(f'{intervention_subcategory}.json')
#      #output_file = os.path.join(output_dir, sanitized_filename)
#      json_file = os.path.join(output_dir, sanitized_filename)
# #     # with open(output_file, 'w') as f:
# #     #     json.dump(json_output, f, indent=2)
# #     #
# #     # print(f"\n\nGenerated JSON file: {sanitized_filename} in the output folder.")
#      with open(json_file, 'w') as f:
#         json.dump(json_output, f, indent = 2)
# 
#      print(f"Generated JSON file: {json_file}")
# 
#      return json_file
# 
# ####### edited indentation
# def generate_output(input_json: str):
#     global biomass_co2_result
#     json_data, df, AOIs, intervention_subcategory = load_json_data(input_json)
# 
#     np.random.seed(1)
#     mc = [GHGcalc(aoi_id, df, 1000, intervention_subcategory, biomass_co2_result) for aoi_id in AOIs]
#     mc = [result for result in mc if not result.empty]  # to filter out empty dataframes
# 
#     if not mc:
#         print("Error: No valid results generated. Check input data.")
#         return None
# 
#     result_list = []
#     for iteration in mc:
#         aoi_id = iteration['AOI'].iloc[0]
#         result_dict = {'AOI': aoi_id, 'Area': iteration['Area'].iloc[0]}
# 
#         for column in iteration.columns:
#             if column not in ['AOI', 'Area']:
#                 result_dict[column] = iteration[column].mean()
#                 result_dict[f"{column}_sd"] = iteration[column].std()
# 
#         result_list.append(result_dict)
# 
#     result_df = pd.DataFrame(result_list)
# 
#     # Calculate totals
#     result_df['N2O_thayr'] = result_df[[col for col in result_df.columns if col.startswith('N2O_') and not col.endswith('_sd')]].sum(axis=1)
#     result_df['sdN2Oha'] = np.sqrt((result_df[[col for col in result_df.columns if col.startswith('N2O_') and col.endswith('_sd')]]**2).sum(axis=1))
#     result_df['CH4_thayr'] = result_df[[col for col in result_df.columns if col.startswith('CH4_') and not col.endswith('_sd')]].sum(axis=1)
#     result_df['sdCH4ha'] = np.sqrt((result_df[[col for col in result_df.columns if col.startswith('CH4_') and col.endswith('_sd')]]**2).sum(axis=1))
#     result_df['CO2_thayr'] = result_df['SOC'] + result_df['totalC']
#     result_df['sdCO2ha'] = np.sqrt(result_df['SOC_sd']**2 + result_df['totalC_sd']**2)
#     result_df['CO2_tyr'] = result_df['CO2_thayr'] * result_df['Area']
#     result_df['sdCO2'] = result_df['sdCO2ha'] * result_df['Area']
#     result_df['N2O_tyr'] = result_df['N2O_thayr'] * result_df['Area']
#     result_df['sdN2O'] = result_df['sdN2Oha'] * result_df['Area']
#     result_df['CH4_tyr'] = result_df['CH4_thayr'] * result_df['Area']
#     result_df['sdCH4'] = result_df['sdCH4ha'] * result_df['Area']
# 
#     # Drop individual columns after totaling
#     columns_to_drop = [col for col in result_df.columns if (col.startswith('N2O_') or col.startswith('CH4_')) and col not in ['N2O_thayr', 'CH4_thayr']] + ['SOC', 'totalC', 'SOC_sd', 'totalC_sd']
#     result_df = result_df.drop(columns=columns_to_drop)
# 
#     # Add total row
#     totals = result_df.sum(numeric_only=True).to_frame().T
#     totals['AOI'] = 'Total'
#     result_df = pd.concat([result_df, totals], ignore_index=True)
# 
#     # Reorder columns
#     cols = ['AOI', 'Area'] + [col for col in result_df.columns if col not in ['AOI', 'Area']]
#     result_df = result_df[cols]
# 
#     # Prepare JSON output
#     json_output = []
#     for _, row in result_df.iterrows():
#         aoi_dict = {
#             'AOI': row['AOI'],
#             'Area': row['Area'],
#             'CO2_thayr': [row['CO2_thayr']] * 20,
#             'sdCO2ha': [row['sdCO2ha']] * 20,
#             'CO2_tyr': [row['CO2_tyr']] * 20,
#             'sdCO2': [row['sdCO2']] * 20,
#             'N2O_thayr': [row['N2O_thayr']] * 20,
#             'sdN2Oha': [row['sdN2Oha']] * 20,
#             'N2O_tyr': [row['N2O_tyr']] * 20,
#             'sdN2O': [row['sdN2O']] * 20,
#             'CH4_thayr': [row['CH4_thayr']] * 20,
#             'sdCH4ha': [row['sdCH4ha']] * 20,
#             'CH4_tyr': [row['CH4_tyr']] * 20,
#             'sdCH4': [row['sdCH4']] * 20
#         }
#         json_output.append(aoi_dict)
# 
#     # Generate output files
#     output_dir = 'output'
#     os.makedirs(output_dir, exist_ok=True)
# 
#     # Write CSV output (optional)
#     csv_file = os.path.join(output_dir, 'result.csv')
#     result_df.to_csv(csv_file, index=False)
#     print(f"Generated CSV file: {csv_file}")
# 
#     # Write JSON output
#     sanitized_filename = sanitize_filename(f'{intervention_subcategory}.json')
#     json_file = os.path.join(output_dir, sanitized_filename)
#     with open(json_file, 'w') as f:
#         json.dump(json_output, f, indent=2)
# 
#     print(f"Generated JSON file: {json_file}")
# 
#     return json_file
# 
# #Testing output
# json_file_path = './ARR_final.json'
# output_file = generate_output(json_file_path)
# print(f"Output file saved at: {output_file}")
# 
# 

## Fixing issues
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

# Example usage
try:
    json_file_path = './ARR_final.json'
    output_file = generate_output(json_file_path)
    if output_file:
        print(f"Output file saved at: {output_file}")
    else:
        print("Failed to generate output file.")
except Exception as e:
    print(f"An error occurred: {str(e)}")
