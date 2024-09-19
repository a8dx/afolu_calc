# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: IFM_final.py
# Author: Barbara Bomfim
# Date Started: 07/14/2024
# Last Edited: 09/19/2024
# Purpose: AFOLU GHG Calculations for Improved Forest Management Interventions (RIL and Extended Rotation)
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data from Dynamic World land use data,  data layers
# User inputs: User inputs depend on sub-intervention category (reduced-impact logging (RIL), Extended Rotation (ER), IFM (stop logging)). 
#   - Initial:
#           -land area: area under business-as-usual land use
#.          -forest_management_type: Natural, Plantation
#.  -Baseline data:
#.          -forest_management_type: Natural, Plantation
#           -forest_type: forest type prior to intervention
#           -D: wood density (t d.m./m3)
#           -VolExt = volume timber over bark extracted before project intervention (m3/ha)
#           -AHA (Annual Harvest Area): in ha - user input. 
#               -- Will attribute a value in the JSON but will indicate this in the code.
#           -C_timber: Proportion of total carbon extracted that resides in each wood product class 
#                 -- (Roundwood: Tropical "0.55", Temperate "0.45", Boreal "0.42", Sawnwood: Tropical 0.58, Temperate 0.48, Boreal 0.44, Woodbase panels: Tropical 0.62, Temperate 0.52, Boreal 0.52)
#           -prop_ox: Carbon emitted due to short-term oxidation of wood products "0.25" - the average across: Sawnwood 0.2; Woodbase panels 0.1; Other industrial roundwood 0.3; Paper and paperboard 0.4)
#           -prop_ox_long: Carbon in additional oxidized fraction "0.947" - the average across: Sawnwood 0.84; Woodbase panels 0.97; Other industrial roundwood 0.99; Paper and paperboard 0.99
#           -RLbau: Number of years in baseline rotation length
#           -GR: Growth rate (m3/ha/yr)
#           -BEF: Biomass Expansion Factor
#           -half_life_wp = Half-life of wood products (years)
#  -Intervention requires:
#           -forest_area: forest area under intervention
#.          -forest_management_type: Natural, Plantation
#           -forest_type: forest type prior to intervention
#           -D: wood density (t/m3)
#           -VolExt = volume timber over bark extracted during project intervention (m3/ha)
#           -AHA (Annual Harvest Area): in ha - user input. 
#               -- Will attribute a value in the JSON but will indicate this in the code.
#           -C_timber: Proportion of total carbon extracted that resides in each wood product class 
#                 -- (Roundwood: Tropical "0.55", Temperate "0.45", Boreal "0.42", Sawnwood: Tropical 0.58, Temperate 0.48, Boreal 0.44, Woodbase panels: Tropical 0.62, Temperate 0.52, Boreal 0.52)
#           -prop_ox: Carbon emitted due to short-term oxidation of wood products "0.25" - the average across: Sawnwood 0.2; Woodbase panels 0.1; Other industrial roundwood 0.3; Paper and paperboard 0.4)
#           -prop_ox_long: Carbon in additional oxidized fraction "0.947" - the average across: Sawnwood 0.84; Woodbase panels 0.97; Other industrial roundwood 0.99; Paper and paperboard 0.99
#           -RLint: Number of years in intervention rotation length
#           -GR: Growth rate (m3/ha/yr)
#           -BEF: Biomass Expansion Factor
#           -half_life_wp = Half-life of wood products (years)

#### Parameters and Paths ####
# Calculations requires the following parameters and their associated uncertainty:

##AGB
#AGBREF:AGB for the climate zone and soil type (tC/ha)

##BGB
#BGBREF:BGB for the climate zone and soil type (tC/ha)

##SOC

##N2O

##CH4

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
import json
import os
import math
from typing import Dict, List, Any # edited this line

import pandas as pd
import numpy as np
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor
import ee
import re

from helpersIFM import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock, sanitize_filename

## Step 0: troubleshooting
# Define function to check for required global variables
def check_global_variables():
    required_globals = ['data', 'polygon', 'average_agbd_dict']
    for var in required_globals:
        if var not in globals():
            raise NameError(f"Required global variable '{var}' is not defined")

### STEP 1 - Estimate mean annual AGB (in Mg/ha) using GEDI dataset ####
# GEE path: https://code.earthengine.google.com/?scriptPath=users%2Fbabomfimf%2FAGB-GEDI-L4A%3Atest-3_AGB_annual_mean
## Anthony's GEE: https://code.earthengine.google.com/44d6aaa5db4764b5b5f3825baf900d04

# Authenticate GEE #
# if this doesn't work authenticate from the command line  by running `earthengine authenticate`
ee.Authenticate()

# Initialize Earth Engine
# Use the GEE project name
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
file_name = 'IFM_final.json'
file_path = os.path.join(current_directory, file_name)
json_file_path = './IFM_final.json'
# Load data from the JSON file
with open(file_path, 'r') as file:
        data = json.load(file)

### STEP 2: Calculate Mean annual Aboveground Carbon (AGC) stock in tCO2e/ha ####
# here we obtain aboveground biomass and convert to aboveground carbon stock
def mean_annual_agc_stock(scenario, log_level='info'):
    global data, polygon, average_agbd_dict
    
    scenario_data = data["scenarios"][scenario][0]
    
    average_agbd_value = average_agbd_dict.get('agbd_mean', 0)
    area_converted = polygon.area().getInfo()/10000  # converting to hectares
    area_converted_yr = area_converted
    
    CF = scenario_data.get("CF", 0)
    average_agbd_cstock = convert_to_c(average_agbd_value, CF)
    average_agbd_tco2e = convert_to_co2e(average_agbd_cstock)
    
    subregion_results = []
    for subregion in scenario_data["aoi_subregions"]:
        subregion_area = subregion["area"]
        subregion_fraction = subregion_area / area_converted_yr if area_converted_yr != 0 else 0
        
        subregion_agbd_tco2e = average_agbd_tco2e * subregion_fraction
        
        subregion_results.append({
            "aoi_id": subregion["aoi_id"],
            "area": subregion_area,
            "carbon_stock": subregion_agbd_tco2e
        })
    
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

### STEP 4: Annual CO2 Calculations for both Baseline and RIL intervention #########
## Equations:
# Logging emissions - baseline and RIL same equations
# Carbon stored in long-term wood products - baseline and RIL same equations
# Carbon emissions from incidental damage - baseline and RIL same equations

## Step 4.1: Define Baseline Conventional Logging Emissions equation
def calculate_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AHA = float(scenario_data.get("AHA", 0)) # this should be user input, as mentioned in User Input commented part of python code
        VolExt = float(scenario_data.get("VolExt", 0))
        D = float(scenario_data.get("D", 0))
        ElE_factor_1 = float(scenario_data.get("ElE_factor_1", 0))
        ElE_factor_2 = float(scenario_data.get("ElE_factor_2", 0))
        
        ELE = (ElE_factor_1 * D) - ElE_factor_2
        TimberTree = (AHA * VolExt * ELE)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "TimberTree": TimberTree,
            "ELE": ELE
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Extracted Log Emissions (ELE): {ELE:.4f} tC/m3 extracted")
            print(f"  Logging emissions : {TimberTree:.2f} tC")
    
    return total_results

## Step 4.2: Define equation to calculate carbon stored in long-term wood products
def calculate_carbon_wood_products(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        C_timber = float(scenario_data.get("C_timber", 0))
        r_w = float(scenario_data.get("r_w", 0)) # this is a standard, fixed value
        prop_ox = 0.25
        #prop_ox = float(scenario_data.get("prop_ox", 0))
        prop_ox_long = 0.947
        #prop_ox_long = float(scenario_data.get("prop_ox_long", 0))
        
        C_ww = C_timber * r_w
        C_slp = (C_timber - C_ww) * prop_ox
        C_addlox = (C_timber - C_ww - C_slp) * prop_ox_long
        
        WoodProducts = C_timber - (C_ww + C_slp + C_addlox)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "WoodProducts": WoodProducts
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Carbon in timber class: {C_timber:.2f} t C")
            print(f"  Carbon in wood waste class: {C_ww:.2f} t C")
            print(f"  Carbon oxidized in less than 5 years: {C_slp:.2f} t C")
            print(f"  Carbon oxidized between 5-100 years: {C_addlox:.2f} t C")
            print(f"  Carbon permanently stored in long term wood products: {WoodProducts:.2f} t C")
    
    return total_results

## Step 4.3: Define equation to calculate Carbon Emissions from Incidental Damage
def calculate_incidental_damage(scenario, log_level='info'):
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
        average_total_tco2e = float(scenario_data.get("average_total_tco2e", 0))
        
        LDF = LDF_factor_1 * average_total_tco2e + LDF_factor_2
        IncidentalDamage = (LDF * VolExt) * 44/12
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "LDF": LDF,
            "IncidentalDamage": IncidentalDamage
        })
    
    return total_results

## Step 4.4: Define equation to calculate Conventional carbon emissions from Logging Infrastructure
def calculate_logging_infrastructure(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AHA = float(scenario_data.get("AHA", 0)) # this should be user input, as mentioned in User Input commented part of python code
        VolExt = float(scenario_data.get("VolExt", 0))
        SkidsFactor = float(scenario_data.get("SkidsFactor", 0))
        RoadsDecksFactor = float(scenario_data.get("RoadsDecksFactor", 0))
        
        Skids = SkidsFactor * AHA * VolExt
        RoadsDeck = RoadsDecksFactor * AHA * VolExt
        LoggingInfrastructure = Skids + RoadsDeck
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "Skids": Skids,
            "RoadsDeck": RoadsDeck,
            "LoggingInfrastructure": LoggingInfrastructure
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Annual Harvest Area (AHA): {AHA:.2f} ha")
            print(f"  Volume Timber Extracted (VolExt): {VolExt:.2f} m3/ha")
            print(f"  Skids emissions: {Skids:.2f} tC")
            print(f"  Roads and Decks Emissions: {RoadsDeck:.2f} tC")
            print(f"  Logging Infrastructure Emissions: {LoggingInfrastructure:.2f} tC")
    
    return total_results

## Step 4.5: Estimate Total Logging Emissions
def calculate_total_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data["aoi_subregions"]
    
    timber_tree_results = calculate_logging_emissions(scenario, log_level)
    wood_products_results = calculate_carbon_wood_products(scenario, log_level)
    incidental_damage_results = calculate_incidental_damage(scenario, log_level)
    logging_infrastructure_results = calculate_logging_infrastructure(scenario, log_level)
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        timber_tree = next(item for item in timber_tree_results if item["aoi_id"] == aoi_id)["TimberTree"]
        wood_products = next(item for item in wood_products_results if item["aoi_id"] == aoi_id)["WoodProducts"]
        incidental_damage = next(item for item in incidental_damage_results if item["aoi_id"] == aoi_id)["IncidentalDamage"]
        logging_infrastructure = next(item for item in logging_infrastructure_results if item["aoi_id"] == aoi_id)["LoggingInfrastructure"]
        
        Logging_Emissions = ((timber_tree - wood_products) + incidental_damage + logging_infrastructure) * (44 / 12)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "TimberTree": timber_tree,
            "WoodProducts": wood_products,
            "IncidentalDamage": incidental_damage,
            "LoggingInfrastructure": logging_infrastructure,
            "Logging_Emissions": Logging_Emissions
        })
    
    return total_results

###### CALCULATE SUB-INTERVENTIONS CO2 IMPACT #######

### Step 5: Calculate RIL variables ####
## Calculate bau and intervention scenarios and difference between both

# def calculate_all(scenario, log_level='info'):
#     result = {}
#     
#     average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
#     result["average_agbd_tco2e"] = average_agbd_tco2e
#     
#     mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
#     result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
#     
#     TimberTree = calculate_logging_emissions(scenario=scenario, log_level=log_level)
#     result["TimberTree"] = TimberTree
#     if log_level == 'debug':
#         print(f"Annual Logging Emissions: {TimberTree} tCO2e/yr")
#     
#     WoodProducts = calculate_carbon_wood_products(scenario=scenario, log_level=log_level)
#     result["WoodProducts"] = WoodProducts
#     if log_level == 'debug':
#         print(f"Carbon stored in wood products: {WoodProducts} tCO2e")
#         
#     IncidentalDamage = calculate_incidental_damage(scenario=scenario, log_level=log_level)
#     result["IncidentalDamage"] = IncidentalDamage
#     if log_level == 'debug':
#         print(f"Logging Emissions due to Incidental Damage: {IncidentalDamage} tCO2e/yr")
# 
#     LoggingInfrastructure = calculate_logging_infrastructure(scenario=scenario, log_level=log_level)
#     result["LoggingInfrastructure"] = LoggingInfrastructure
#     if log_level == 'debug':
#         print(f"Logging Emissions due to Incidental Damage: {LoggingInfrastructure} tCO2e/yr")
# 
#     Logging_Emissions = calculate_total_logging_emissions(scenario=scenario, log_level=log_level)
#     result["Logging_Emissions"] = Logging_Emissions
#     if log_level == 'debug':
#         print(f"Logging Emissions due to Incidental Damage: {Logging_Emissions} tCO2e/yr")
# 
#     return result
def calculate_all(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    logging_emissions = calculate_total_logging_emissions(scenario=scenario, log_level=log_level)
    result["Logging_Emissions"] = logging_emissions
    
    return result

# # biomass = {}
# # biomass["business_as_usual"] = calculate_all(scenario="business_as_usual")
# # biomass["intervention"] = calculate_all(scenario="intervention")
# # 
# # biomass_co2_result = {}
# # 
# # for inter_item, bau_item in zip(biomass["intervention"]["Logging_Emissions"], 
# #                                 biomass["business_as_usual"]["Logging_Emissions"]):
# #     if inter_item['aoi_id'] != bau_item['aoi_id']:
# #         raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
# #     
# #     aoi_id = inter_item['aoi_id']
# #     difference = inter_item['Logging_Emissions'] - bau_item['Logging_Emissions']
# #     
# #     biomass_co2_result[aoi_id] = difference# inter_result = calculate_all(scenario="intervention")
# # 
# # print(biomass_co2_result)
# # 
# # biomass_co2_result_error_positive = 10  # Assuming this is a percentage
# # 
# # biomass_co2_sd = {}
# # 
# # for aoi_id, difference in biomass_co2_result.items():
# #     # Calculate SD for each aoi_id
# #     sd = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96
# #     biomass_co2_sd[aoi_id] = sd

biomass = {
    "business_as_usual": calculate_all(scenario="business_as_usual"),
    "intervention": calculate_all(scenario="intervention")
}

biomass_co2_result = {}
biomass_co2_sd = {}
biomass_co2_result_error_positive = 10  # this is a percentage

def get_annual_change(scenario_result):
    if 'Logging_Emissions' in scenario_result:
        return scenario_result['Logging_Emissions']
    else:
        print(f"Warning: 'Logging_Emissions' not found in {scenario_result.get('scenario', 'unknown')} scenario.")
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
    
    if 'Logging_Emissions' not in inter_item or 'Logging_Emissions' not in bau_item:
        print(f"Warning: Missing Logging_Emissions for aoi_id {aoi_id}. Skipping.")
        continue

    difference = inter_item['Logging_Emissions'] - bau_item['Logging_Emissions']
    
    biomass_co2_result[aoi_id] = difference
    biomass_co2_sd[aoi_id] = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96

print("Biomass CO2 Result:", biomass_co2_result)
print("Biomass CO2 Standard Deviation:", biomass_co2_sd)

## To Mathematica Team: Check why output is like this:
# print("Biomass CO2 Result:", biomass_co2_result)
# Biomass CO2 Result: {'tropicalmoist_acrisols': -7030.206333333328, 'tropicalmoist_ferralsols': -7030.206333333328}
# >>> print("Biomass CO2 Standard Deviation:", biomass_co2_sd)
# Biomass CO2 Standard Deviation: {'tropicalmoist_acrisols': 358.6839965986392, 'tropicalmoist_ferralsols': 358.6839965986392}


####### EXTENDED ROTATION INTERVENTION ##########
## Defining equations for use if Extended Rotation is the selected sub-intervention

## Step 1-ER: Define equation to calculate baseline total tree carbon stock
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
  
## Step 2-ER: Define equation to calculate wood product carbon stock
def calculate_wood_products_c_stock(scenario, t=None, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        average_agbd_tco2e = subregion["carbon_stock"]
        BEF = float(scenario_data.get("BEF", 0)) # Biomass expansion factor for wood removal
        ConvEff = 0.5 # to check: milling conversion efficiency (ConvEff)
        
        WPBL0 = average_agbd_tco2e / BEF * ConvEff
        
        result = {
            "aoi_id": aoi_id,
            "area": area,
            "WPBL0": WPBL0
        }
        
        if t is not None:
            half_life_wp = float(scenario_data.get("half_life_wp", 0))
            WPBLt = math.exp(-((math.log(0.5) / half_life_wp) * t)) * WPBL0
            result["WPBLt"] = WPBLt
            
            if log_level == 'debug':
                print(f"Subregion {aoi_id}:")
                print(f"  Wood products at time t baseline: {WPBLt:.2f} t CO2e")
        
        total_results.append(result)
        
        if log_level == 'debug' and t is None:
            print(f"Subregion {aoi_id}:")
            print(f"  Wood products at time zero baseline: {WPBL0:.2f} t CO2e")
    
    return total_results

### Step 3-ER: Define equation to calculate emissions from Conventional Logging as following:
## Emissions (t CO2e) = Harvest area * (result from Step 1 + result from Step 2)
def calculate_conventional_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    
    # Get results from Step 1-ER:
    subregion_results = mean_annual_agc_stock(scenario, log_level)  # Assuming this function exists
    total_c_stock_results = mean_annual_tot_c_stock(subregion_results, scenario, log_level)
    
    # Get results from Step 2-ER:
    wood_products_results = calculate_wood_products_c_stock(scenario, log_level=log_level)
    
    emissions_results = []
    
    for c_stock, wood_product in zip(total_c_stock_results, wood_products_results):
        aoi_id = c_stock['aoi_id']
        area = c_stock['area']
        
        # Get the harvest area from scenario data or use a default value
        harvest_area = float(scenario_data.get("harvest_area", area))
        
        # Calculate emissions
        total_c_stock = c_stock['total_carbon_stock']
        wood_products_c_stock = wood_product['WPBL0']
        
        emissions = harvest_area * (total_c_stock - wood_products_c_stock)
        
        emissions_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "harvest_area": harvest_area,
            "total_carbon_stock": total_c_stock,
            "wood_products_carbon_stock": wood_products_c_stock,
            "emissions": emissions
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Harvest Area: {harvest_area:.2f} ha")
            print(f"  Total Carbon Stock: {total_c_stock:.2f} tCO2e/ha")
            print(f"  Wood Products Carbon Stock: {wood_products_c_stock:.2f} tCO2e/ha")
            print(f"  Emissions from Conventional Logging: {emissions:.2f} tCO2e")
    
    return emissions_results

##### Step 4-ER: Define equations to calculate ER GHG impact ###

## Step 4.1-ER: Define equation to calculate Biomass Accumulation due to ER
def calculate_extended_rotation_live_biomass(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        average_agbd_tco2e = subregion["carbon_stock"]
        average_bgbd_tco2e = convert_to_bgb(average_agbd_tco2e, ratio)
        average_total_tco2e = average_agbd_tco2e + average_bgbd_tco2e
        
        RLbau = float(scenario_data.get("RLbau", 0)) # to check: rotation length bau
        RLint = float(scenario_data.get("RLint", 0)) #to check: rotation length int
        GR = float(scenario_data.get("GR", 0))
        RS = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))
        
        ERY = RLint - RLbau
        ERB = (average_agbd_tco2e + (ERY * GR)) * (1 + RS)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "ERB": ERB
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"Extended rotation live biomass: {ERB:.2f} t CO2e")
    
    return total_results

## Step 4.2: Define equation to calculate difference in years between RLbau and RLint (ER)
def calculate_ery(json_file_path):
    try:
        # Load JSON data
        with open(json_file_path, 'r') as file:
            data = json.load(file)
    except json.JSONDecodeError as e:
        print(f"JSON decoding error: {str(e)}")
        print("Please check your JSON file for formatting errors.")
        return None
    except FileNotFoundError:
        print(f"File not found: {json_file_path}")
        return None
    
    # Extract RL values from scenarios
    scenarios = data.get('scenarios', {})
    bau_rl = None
    int_rl = None
    
    if 'business_as_usual' in scenarios and scenarios['business_as_usual']:
        bau_rl = scenarios['business_as_usual'][0].get('RL')
    if 'intervention' in scenarios and scenarios['intervention']:
        int_rl = scenarios['intervention'][0].get('RL')
    
    # Check if both RL values are present
    if bau_rl is None:
        print("Missing RL value for 'business_as_usual' scenario")
        return None
    if int_rl is None:
        print("Missing RL value for 'intervention' scenario")
        return None
    
    # Ensure RL values are numeric
    try:
        bau_rl = float(bau_rl)
        int_rl = float(int_rl)
    except ValueError:
        print("RL values must be numeric")
        return None
    
    # Calculate ERY
    ery = bau_rl - int_rl
    
    return {
        'business_as_usual_RL': bau_rl,
        'intervention_RL': int_rl,
        'ERY': ery
    }

# testing
json_file_path = 'IFM_final.json'
result = calculate_ery(json_file_path)
if result:
    print(f"Business-as-usual RL: {result['business_as_usual_RL']}")
    print(f"Intervention RL: {result['intervention_RL']}")
    print(f"Emissions Reduction per Year (ERY): {result['ERY']}")
else:
    print("Failed to calculate ERY. Please check the error messages above.")


### Step 4.3-ER: Define equation to calculate carbon in wood products due to ER
def calculate_wood_products_extended_rotation(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        average_agbd_tco2e = subregion["carbon_stock"]
        RLbau = float(scenario_data.get("RLbau", 0)) # to check
        RLint = float(scenario_data.get("RLint", 0)) # to check
        GR = float(scenario_data.get("GR", 0)) #growth rate
        BEF = float(scenario_data.get("BEF", 0)) #Biomass Expansion Factor
        
        ERY = RLint - RLbau
        WPER0 = (average_agbd_tco2e + (ERY * GR)) / BEF
        
        result = {
            "aoi_id": aoi_id,
            "area": area,
            "WPER0": WPER0
        }
        
        if t is not None:
            half_life_wp = float(scenario_data.get("half_life_wp", 0))
            WPERt = math.exp(-((math.log(0.5) / half_life_wp) * t)) * WPER0
            result["WPERt"] = WPERt
            
            if log_level == 'debug':
                print(f"Subregion {aoi_id}:")
                print(f"  Wood products at time t extended rotation: {WPERt:.2f} t CO2e")
        
        total_results.append(result)
        
        if log_level == 'debug' and t is None:
            print(f"Subregion {aoi_id}:")
            print(f"  Wood products at time zero extended rotation: {WPER0:.2f} t CO2e")
    
    return total_results

#Define equation to calculate the live biomass benefit of extended rotation
def calculate_benefit_extended_rotation_live_biomass(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    ratio = scenario_data["ratio_below_ground_biomass_to_above_ground_biomass"]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        RLbau = float(scenario_data.get("RLbau", 0))
        RLint = float(scenario_data.get("RLint", 0))
        GR = float(scenario_data.get("GR", 0))
        RS = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))
        
        average_agbd_tco2e = subregion["carbon_stock"]
        average_bgbd_tco2e = convert_to_bgb(average_agbd_tco2e, ratio)
        average_total_tco2e = average_agbd_tco2e + average_bgbd_tco2e
        
        ERY = RLint - RLbau
        ERB = (average_agbd_tco2e + (ERY * GR)) * (1 + RS)
        
        benefit_live_biomass_ER = ((ERB * area) - (average_total_tco2e * area)) * (44 / 12)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "benefit_live_biomass_ER": benefit_live_biomass_ER
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Total benefit long term live biomass extended rotation: {benefit_live_biomass_ER:.2f} t CO2e")
    
    return total_results

#Define equation to calculate the wood products benefit of extended rotation
def calculate_benefit_extended_rotation_wood_products(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        benefit_wood_products_ER = ((WPERt * area) - (WPER0 * area)) * (44 / 12)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "benefit_wood_products_ER": benefit_wood_products_ER
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Benefit wood products extended rotation: {benefit_wood_products_ER:.2f} t CO2e")
    
    return total_results

def calculate_total_benefits(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        benefit_live_biomass_ER = float(scenario_data.get("benefit_live_biomass_ER", 0))
        benefit_wood_products_ER = float(scenario_data.get("benefit_wood_products_ER", 0))
        
        total_benefits = benefit_live_biomass_ER + benefit_wood_products_ER
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "total_benefits": total_benefits
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Total benefits: {total_benefits:.2f} t CO2e")
    
    return total_results

#####Calculate all IFM - Extended Rotation variables ####
def calculate_all_er(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    BL = calculate_baseline_biomass_carbon_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["BL"] = BL
    
    WPB = calculate_wood_products_baseline(scenario=scenario, log_level=log_level)
    result["WPB"] = WPBL0
    if log_level == 'debug':
        print(f"Wood products baseline: {WPB} tCO2e")
        
    Long_live_biomass = calculate_long_term_live_biomass(scenario=scenario, log_level=log_level)
    result["Long_live_biomass"] = Long_live_biomass
    if log_level == 'debug':
        print(f"Long term live tree biomass: {Long_live_biomass} tCO2e")
    
    WoodProducts_0 = calculate_wood_products_extended_rotation_zero(scenario=scenario, log_level=log_level)
    result["WoodProducts_0"] = WoodProducts
    if log_level == 'debug':
        print(f"Carbon stored in wood products ER time ze3ro: {WoodProducts_0} tCO2e")
    
    WoodProducts = calculate_wood_products_extended_rotation_time_t(scenario=scenario, log_level=log_level)
    result["WoodProducts"] = WoodProducts
    if log_level == 'debug':
        print(f"Carbon stored in wood products ER: {WoodProducts} tCO2e")
        
    Biomass_wood_benefits = calculate_benefit_extended_rotation_live_biomass(scenario=scenario, log_level=log_level)
    result["Biomass_wood_benefits"] = Biomass_wood_benefits
    if log_level == 'debug':
        print(f"Live biomass benefits of ER: {Biomass_wood_benefits} tCO2e/yr")

    ER_wood_benefits = calculate_benefit_extended_rotation_wood_products(scenario=scenario, log_level=log_level)
    result["ER_wood_benefits"] = ER_wood_benefits
    if log_level == 'debug':
        print(f"Wood products benefits of ER: {ER_wood_benefits} tCO2e/yr")

    ER_total_benefits = calculate_total_benefits(scenario=scenario, log_level=log_level)
    result["ER_total_benefits"] = ER_total_benefits
    if log_level == 'debug':
        print(f"Total benefits of ER: {ER_total_benefits} tCO2e/yr")

    return result

biomass_er = {}
biomass_er["business_as_usual"] = calculate_all_er(scenario="business_as_usual")
biomass_er["intervention"] = calculate_all_er(scenario="intervention")

biomass_er_co2_result = {}

for inter_item, bau_item in zip(biomass_er["intervention"]["ER_total_benefits"], 
                                biomass_er["business_as_usual"]["ER_total_benefits"]):
    if inter_item['aoi_id'] != bau_item['aoi_id']:
        raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
    
    aoi_id = inter_item['aoi_id']
    difference = inter_item['ER_total_benefits'] - bau_item['ER_total_benefits']
    
    biomass_er_co2_result[aoi_id] = difference inter_result = calculate_all(scenario="intervention")

print(biomass_er_co2_result)

biomass_er_co2_result_error_positive = 10  # Assuming this is a percentage

biomass_er_co2_sd = {}

for aoi_id, difference in biomass_er_co2_result.items():
    # Calculate SD for each aoi_id
    sd = (abs(difference) * biomass_er_co2_result_error_positive / 100) / 1.96
    biomass_er_co2_sd[aoi_id] = sd


######## Getting JSON data ready to be calculated ######
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
    
    # Remove other uncertainty columns
    df = df.drop(columns=[col for col in df.columns if 'error' in col])
    
    return json_data, df, AOIs, interv_sub

###############################################################

def average_results_sims(input_dict):
    results_unc = {}
    for key, value in input_dict.items():
        if isinstance(value, list):
            if value:  # Check if the list is not empty
                average = sum(value) / len(value)
                results_unc[key] = average
            else:
                results_unc[key] = 0  # or None, for empty lists
        else:
            results_unc[key] = value  # Keep non-list values as they are
    return results_unc

def GHGcalc(aoi_id, df, nx, intervention_subcategory, biomass_co2_result):
    
    temp_bau = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "business-as-usual")]
    temp_int = df[(df['aoi_id'] == aoi_id) & (df['scenario'] == "intervention")]
    # Create a structured array that can hold both strings and floats
    results_unc = {}
    
    results_unc['AOI'] = aoi_id
    results_unc['Area'] = temp_bau['area'].values[0]
    #results_unc['SOC'] = [] ## SOC not included
    results_unc['totalC'] = []
    results_unc[f'N2O_{aoi_id}'] = []
    # Note: CH4 calculation not included
    # results_unc[f'CH4_{aoi_id}'] = []
    
    for m in range(nx):
        # SOC Stock Change not included
        # here we add the biomass specific to the aoi_id
        results_unc["totalC"].append(biomass_co2_result[aoi_id]) 
       
        # results_unc[f"N2O_{aoi_id}"].append(N2O)
      

    # column_names = ['AOI', 'Area', 'Rep', 'SOC', 'totalC', f'N2O_y{aoi_id}', f'CH4_y{aoi_id}']
    print("UNC", results_unc)
    results_unc_df = pd.DataFrame(results_unc)
    return results_unc_df

def generate_output(input_json):
    global biomass_co2_result
    json_data, df, AOIs, intervention_subcategory = load_json_data(input_json)
    
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
    # Run code below to roll them up for totals as well instead of aoi_id
    # Index(['AOI', 'Area', 'Rep', 'SOC', 'SOC_sd', 'totalC', 'totalC_sd',
    #    'N2O_tropicalmoist_acrisols', 'N2O_tropicalmoist_acrisols_sd',
    #    'CH4_tropicalmoist_acrisols', 'CH4_tropicalmoist_acrisols_sd',
    #    'N2O_tropicalmoist_ferralsols', 'N2O_tropicalmoist_ferralsols_sd',
    #    'CH4_tropicalmoist_ferralsols', 'CH4_tropicalmoist_ferralsols_sd'],
    #   dtype='object')
    
    # Summing specific columns
    #result_df['N2O_thayr'] = result_df['N2O_tropicalmoist_acrisols'] + result_df['N2O_tropicalmoist_ferralsols']
    #result_df['sdN2Oha'] = result_df['N2O_tropicalmoist_acrisols_sd'] + result_df['N2O_tropicalmoist_ferralsols_sd']
    # result_df['CH4_thayr'] = result_df['CH4_tropicalmoist_acrisols'] + result_df['CH4_tropicalmoist_ferralsols']
    # result_df['sdCH4ha'] = result_df['CH4_tropicalmoist_acrisols_sd'] + result_df['CH4_tropicalmoist_ferralsols_sd']
    result_df['CO2_thayr'] = result_df['SOC'] + result_df['totalC']
    result_df['sdCO2ha'] = result_df['SOC_sd'] + result_df['totalC_sd']
    result_df['CO2_tyr'] = result_df['CO2_thayr'] * result_df['Area']
    result_df['sdCO2'] = result_df['sdCO2ha'] * result_df['Area']
    #result_df['N2O_tyr'] = result_df['N2O_thayr'] * result_df['Area']
    #result_df['sdN2O'] = result_df['sdN2Oha'] * result_df['Area']
    # result_df['CH4_tyr'] = result_df['CH4_thayr'] * result_df['Area']
    # result_df['sdCH4'] = result_df['sdCH4ha'] * result_df['Area']

    # Drop the individual columns after summing
    result_df = result_df.drop(columns=[
        'totalC', 'totalC_sd',
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
                   ['sdCO2ha'] + [f'sdN2Oha_y{i}' for i in range(1, 21)] + \
                   ['CO2_tyr'] + [f'N2O_tyr_y{i}' for i in range(1, 21)] + \
                   ['sdCO2'] + [f'sdN2O_y{i}' for i in range(1, 21)]
                # add: CH4 calculation to be added
                   # [f'CH4_thayr_y{i}' for i in range(1, 21)] + \
                   # [f'sdCH4ha_y{i}' for i in range(1, 21)] + \
                   # [f'CH4_tyr_y{i}' for i in range(1, 21)] + \
                   # [f'sdCH4_y{i}' for i in range(1, 21)]

    resdf = pd.DataFrame(result_df, columns=column_names)

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
                #TODO: CH4 calculation needed
            # 'CH4_thayr': resdf.iloc[r][[f'CH4_thayr_y{i}' for i in range(1, 21)]].tolist(),
            # 'sdCH4ha': resdf.iloc[r][[f'sdCH4ha_y{i}' for i in range(1, 21)]].tolist(),
            # 'CH4_tyr': resdf.iloc[r][[f'CH4_tyr_y{i}' for i in range(1, 21)]].tolist(),
            # 'sdCH4': resdf.iloc[r][[f'sdCH4_y{i}' for i in range(1, 21)]].tolist()
        }
        restib.append(aoi_dict)

    # Write the DataFrame to a CSV file
    result_df.to_csv('result.csv', index=False)
    # 
    json_output = json.dumps(result_df.to_dict(orient='records'), indent=2)
    # 
    # # # Ensure the output directory exists
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    # 
    # # Write the JSON output
    sanitized_filename = 'output_file.json'  # Replace with actual filename sanitization if needed
    output_file = os.path.join(output_dir, sanitized_filename)
    with open(output_file, 'w') as f:
        json.dump(json_data, f, indent=4)

    sanitized = sanitize_filename(f'{intervention_subcategory}.json')
    # 
    # # # Write the JSON output
    output_file = os.path.join(output_dir, sanitized)
    with open(output_file, 'w') as f:
       f.write(json_output)

    print(f"\n\nGenerated file: {sanitized_filename} in the output folder.")

generate_output(json_file_path)
