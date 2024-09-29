# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: IFM_final.py
# Author: Barbara Bomfim
# Date Started: 07/14/2024
# Last Edited: 09/29/2024
# Purpose: AFOLU GHG Calculations for Improved Forest Management Interventions (Reduced-Impact Logging, Extended Rotation, Stop Logging, Fuelwood)
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
#           -AHA_perc (Annual Harvest Area percentage): the percentage of the total area that is actually harvestd. In Json, 0.5.
#           -RL: Number of years in baseline rotation length
#  -Intervention requires:
#           -forest_area: forest area under intervention
#.          -forest_management_type: Natural, Plantation
#           -forest_type: forest type prior to intervention
#           -D: wood density (t/m3)
#           -VolExt = volume timber over bark extracted during project intervention (m3/ha)
#           -AHA_perc (Annual Harvest Area percentage): the percentage of the total area that is actually harvestd. In Json, 0.5.
#           -RL: Number of years in intervention rotation length

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
    
    CF = scenario_data.get("CF", 0)
    average_agbd_cstock = convert_to_c(average_agbd_value, CF)
    average_agbd_tco2e = convert_to_co2e(average_agbd_cstock)
    
    subregion_results = []
    for subregion in scenario_data["aoi_subregions"]:
        subregion_area = subregion["area"]
        subregion_fraction = subregion_area / area_converted if area_converted != 0 else 0
        
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


################ STEP 4 INTERVENTION CALCULATIONS #############################


####### FUELWOOD CONSUMPTION SUB-INTERVENTION #######

### Step 4.1-FW: 
## Define Lfuelwood equation to estimate carbon loss due to fuelwood consumption ###
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
        R = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
        D = scenario_data.get("D")
        CF = scenario_data.get("CF")
        
        if None in [FGtrees, FGpart, BCEFr, R, D, CF]:
            if log_level == 'debug':
                print(f"Missing data for calculating Lfuelwood in subregion {aoi_id}. Please check your data for scenario {scenario}.")
            continue
        
        lfuelwood = ((FGtrees * BCEFr * (1 + R)) + FGpart * D) * CF * (44 / 12)
        
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

### Step 4.2-FW: Calculate Difference between baseline and bau fuelwood consumption

# Function to calculate fuelwood emissions for baseline and intervention scenarios
def calculate_all(scenario, log_level='info'):
    result = {}
    
    try:
        check_global_variables()
        
        steps = [
            ('average_agbd_tco2e', lambda: mean_annual_agc_stock(scenario=scenario, log_level=log_level)),
            ('mean_annual_tot_c_stock', lambda: mean_annual_tot_c_stock(result['average_agbd_tco2e'], scenario=scenario, log_level=log_level)),
            ('lfuelwood', lambda: calculate_lfuelwood(scenario=scenario, log_level=log_level)),
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

## Define calculations to calculate difference between the impact of BAU and intervention scenarios
biomass = {
    "business_as_usual": calculate_all(scenario="business_as_usual"),
    "intervention": calculate_all(scenario="intervention")
}

biomass_co2_result = {}
biomass_co2_sd = {}
biomass_co2_result_error_positive = 10  # this is a percentage

def get_annual_change(scenario_result):
    if 'lfuelwood' in scenario_result:
        return scenario_result['lfuelwood']
    else:
        print(f"Warning: 'lfuelwood' not found in {scenario_result.get('scenario', 'unknown')} scenario.")
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

    difference_fuelwood = inter_item['annual_change'] - bau_item['annual_change']
    
    biomass_co2_result[aoi_id] = difference_fuelwood
    biomass_co2_sd[aoi_id] = (abs(difference_fuelwood) * biomass_co2_result_error_positive / 100) / 1.96

print("Biomass CO2 Result Fuelwood:", biomass_co2_result)
print("Biomass CO2 Standard Deviation Fuelwood:", biomass_co2_sd)



######### REDUCED-IMPACT LOGGING (RIL) SUB-INTERVENTION ##########


### STEP 4-RIL: Annual CO2 Calculations for both Baseline and RIL intervention #########

## Equations:
# Logging emissions - baseline and RIL same equations
# Incidental damage emissions - baseline and RIL same equations

## Step 4.1-RIL: Define Baseline Conventional Logging Emissions equation
def calculate_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AHA_perc = float(scenario_data.get("AHA_perc", 0))
        VolExt = float(scenario_data.get("VolExt", 0))
        D = float(scenario_data.get("D", 0))
        ElE_factor_1 = float(scenario_data.get("ElE_factor_1", 0))
        ElE_factor_2 = float(scenario_data.get("ElE_factor_2", 0))
        
        AHA = area * AHA_perc
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


## Step 4.2-RIL: Define equation to calculate Carbon Emissions from Incidental Damage
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

## Step 4.3-RIL: Define equation to calculate Conventional carbon emissions from Logging Infrastructure
def calculate_logging_infrastructure(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AHA_perc = float(scenario_data.get("AHA_perc", 0)) # this should be user input, as mentioned in User Input commented part of python code
        VolExt = float(scenario_data.get("VolExt", 0))
        SkidsFactor = float(scenario_data.get("SkidsFactor", 0))
        RoadsDecksFactor = float(scenario_data.get("RoadsDecksFactor", 0))
        
        AHA = AHA_perc * area ## edited this to accomodate suggetion
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

## Step 4.4-RIL: Estimate Total Logging Emissions
def calculate_total_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data["aoi_subregions"]
    
    timber_tree_results = calculate_logging_emissions(scenario, log_level)
    incidental_damage_results = calculate_incidental_damage(scenario, log_level)
    logging_infrastructure_results = calculate_logging_infrastructure(scenario, log_level)

    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        timber_tree = next(item for item in timber_tree_results if item["aoi_id"] == aoi_id)["TimberTree"]
        incidental_damage = next(item for item in incidental_damage_results if item["aoi_id"] == aoi_id)["IncidentalDamage"]
        logging_infrastructure = next(item for item in logging_infrastructure_results if item["aoi_id"] == aoi_id)["LoggingInfrastructure"]

        Logging_Emissions = (timber_tree + incidental_damage + logging_infrastructure) * (44 / 12)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "TimberTree": timber_tree,
            "IncidentalDamage": incidental_damage,
            "LoggingInfrastructure": logging_infrastructure,
            "Logging_Emissions": Logging_Emissions
        })
    
    return total_results

###### CALCULATE RIL SUB-INTERVENTION CO2 IMPACT #######

### Step 5-RIL: Calculate RIL variables ####
## Calculate bau and intervention scenarios and difference between both
def calculate_all(scenario, log_level='info'):
    result = {}
    
    try:
        check_global_variables()
        
        steps = [
            ('average_agbd_tco2e', lambda: mean_annual_agc_stock(scenario=scenario, log_level=log_level)),
            ('mean_annual_tot_c_stock', lambda: mean_annual_tot_c_stock(result['average_agbd_tco2e'], scenario=scenario, log_level=log_level)),
            ('logging_emissions', lambda: calculate_total_logging_emissions(scenario=scenario, years=20, average_agbd_tco2e=result['average_agbd_tco2e'], log_level=log_level)),
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

#### Step 6-RIL: Calculate bau and intervention scenarios and difference between both

biomass = {
    "business_as_usual": calculate_all(scenario="business_as_usual"),
    "intervention": calculate_all(scenario="intervention")
}

biomass_co2_result = {}
biomass_co2_sd = {}
biomass_co2_result_error_positive = 10  # this is a percentage

def get_annual_change(scenario_result):
    if 'logging_emissions' in scenario_result:
        return scenario_result['logging_emissions']
    else:
        print(f"Warning: 'logging_emissions' not found in {scenario_result.get('scenario', 'unknown')} scenario.")
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
        print(f"Warning: Missing logging_emissions for aoi_id {aoi_id}. Skipping.")
        continue

    difference_RIL = inter_item['annual_change'] - bau_item['annual_change']
    
    biomass_co2_result[aoi_id] = difference_RIL
    biomass_co2_sd[aoi_id] = (abs(difference_RIL) * biomass_co2_result_error_positive / 100) / 1.96

print("Biomass CO2 Result RIL:", biomass_co2_result)
print("Biomass CO2 Standard Deviation RIL:", biomass_co2_sd)

## To Mathematica Team: Check why output is like this:

# print("Biomass CO2 Result:", biomass_co2_result)
# Biomass CO2 Result: {'tropicalmoist_acrisols': -7030.206333333328, 'tropicalmoist_ferralsols': -7030.206333333328}
# >>> print("Biomass CO2 Standard Deviation:", biomass_co2_sd)
# Biomass CO2 Standard Deviation: {'tropicalmoist_acrisols': 358.6839965986392, 'tropicalmoist_ferralsols': 358.6839965986392}


####### EXTENDED ROTATION SUB-INTERVENTION ##########
## Includes equations to estimate ER benefits related to biomass changes

## Step 4.1-ER: Define equation to calculate baseline total tree carbon stock
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
  

### Step 4.5-ER: Define equation to calculate emissions from Conventional Logging
## Emissions (t CO2e) = Harvest area * (result from Step 4.1 + result from Step 4.2)
def calculate_conventional_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = subregion["forest_area"]
        
        AHA_perc = float(scenario_data.get("AHA_perc", 0)) # this should be user input, as mentioned in User Input commented part of python code

        AHA = forest_area * AHA_perc ## edited this here to accomodate suggestion

    # Get results from Step 1-ER:
    subregion_results = mean_annual_agc_stock(scenario, log_level)
    total_c_stock_results = mean_annual_tot_c_stock(subregion_results, scenario, log_level)
    
    # # Get results from Step 2-ER:
    # wood_products_results = calculate_wood_products_c_stock(scenario, log_level=log_level)
    
    emissions_results = []
    
    for c_stock in zip(total_c_stock_results):
        aoi_id = c_stock['aoi_id']
        forest_area = c_stock['area']
        
        # Get the harvest area from scenario data or use a default value
        AHA_perc = float(scenario_data.get("AHA_perc", area))
        
        # Calculate emissions
        total_c_stock = c_stock['total_carbon_stock']
        
        emissions = (AHA_perc * forest_area) * total_c_stock
        
        emissions_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "harvest_area": harvest_area,
            "total_carbon_stock": total_c_stock,
            "emissions": emissions
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Harvest Area: {harvest_area:.2f} ha")
            print(f"  Total Carbon Stock: {total_c_stock:.2f} tCO2e/ha")
            print(f"  Emissions from Conventional Logging: {emissions:.2f} tCO2e")
    
    return emissions_results

##### Step 5-ER: Define equations to calculate ER CO2 impact ###

## Step 5.1: Define equation to calculate difference in years between BAU and ER rotation length
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
    ery = int_rl - bau_rl
    
    return {
        'business_as_usual_RL': bau_rl,
        'intervention_RL': int_rl,
        'ERY': ery
    }

### Step 5.2-ER: Define equation to calculate Biomass Accumulation due to ER

# ## Step 5.2.1: Define equation to obtain Gtotal (gain in carbon per hectare per year)
# ## forest_type: broadleaf slow-growing, broadleaf fast-growing, conifer, eucalyptus
# ## age_class: primary (primary = natural), secondary > 20 years, secondary < 20 years (secondary = plantation)

def calculate_average_annual_biomass_growth(scenario, log_level='info'):
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = float(subregion["area"])
        forest_type = subregion.get("forest_type", "").lower()
        
        # Retrieve common parameters from scenario data
        R = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))
        CF = float(scenario_data.get("CF", 0))
        
        # Calculate Gtotal based on forest type
        if forest_type in ["broadleaf slow-growing", "broadleaf fast-growing"]:
            Gw = float(scenario_data.get("Gw", 0))
            Gtotal = Gw * (1 + R) * CF * 44/12
        elif forest_type in ["conifer", "eucalyptus"]:
            Iv = float(scenario_data.get("Iv", 0))
            BCEFi = float(scenario_data.get("BCEFi", 0))
            Gtotal = (Iv * BCEFi) * (1 + R) * CF * 44/12
        else:
            raise ValueError(f"Unknown forest type: {forest_type}")
        
        subregion_result = {
            "aoi_id": aoi_id,
            "area": area,
            "forest_type": forest_type,
            "Gtotal": Gtotal
        }
        
        total_results.append(subregion_result)
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area: {area:.2f} ha")
            print(f"  Forest Type: {forest_type}")
            print(f"  R: {R:.4f}")
            print(f"  CF: {CF:.4f}")
            if forest_type in ["broadleaf slow-growing", "broadleaf fast-growing"]:
                print(f"  Gw: {Gw:.4f} tonnes d.m. ha-1 yr-1")
            else:
                print(f"  Iv: {Iv:.4f} m3 ha-1 yr-1")
                print(f"  BCEFi: {BCEFi:.4f}")
            print(f"  Gtotal: {Gtotal:.4f} tonnes CO2e ha-1 yr-1")
    
    if log_level == 'debug':
        total_area = sum(result['area'] for result in total_results)
        total_Gtotal = sum(result['area'] * result['Gtotal'] for result in total_results)
        average_Gtotal = total_Gtotal / total_area if total_area > 0 else 0
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Average Gtotal across all subregions: {average_Gtotal:.4f} tonnes CO2e ha-1 yr-1")
    
    return total_results

## Step 5.2.2-ER: Define equation to calculate gain in biomass carbon stock due to extended rotation
def calculate_extended_rotation_live_biomass(scenario, log_level='info'):
    """
    Calculate the gain in carbon stock due to extended rotation for each subregion.
    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    list: Updated total_results with Extended Rotation Biomass (ERB) in tCO2e for each subregion.
    """
    global data, total_results  # Assuming total_results is saved globally from step 4.2.1
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    # Calculate ERY
    ery_result = calculate_ery(scenario)
    if ery_result is None:
        print("Failed to calculate ERY. Cannot proceed with extended rotation calculation.")
        return None
    
    ERY = ery_result['ERY']
    
    for result in total_results:
        aoi_id = result["aoi_id"]
        area = float(result["area"])
        Gtotal = result["Gtotal"]
        
        subregion = next((s for s in aoi_subregions if s["aoi_id"] == aoi_id), None)
        if subregion is None:
            print(f"Warning: No subregion data found for AOI {aoi_id}. Skipping this subregion.")
            continue
        
        average_agbd_tco2e = float(subregion.get("carbon_stock", 0))
        R = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))
        
        ERB = (average_agbd_tco2e + (ERY * Gtotal)) * (1 + R)
        
        # Add ERB to the existing result dictionary
        result["ERB"] = ERB
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Area: {area:.2f} ha")
            print(f"  Average AGBD: {average_agbd_tco2e:.2f} tCO2e/ha")
            print(f"  ERY: {ERY:.2f} years")
            print(f"  Gtotal: {Gtotal:.4f} tCO2e/ha/year")
            print(f"  Root-to-Shoot Ratio (R): {R:.2f}")
            print(f"  Extended Rotation Biomass (ERB): {ERB:.2f} tCO2e")
    
    if log_level == 'debug':
        total_area = sum(result['area'] for result in total_results)
        total_ERB = sum(result.get('ERB', 0) for result in total_results)
        average_ERB = total_ERB / total_area if total_area > 0 else 0
        print(f"\nTotal area across all subregions: {total_area:.2f} ha")
        print(f"Total Extended Rotation Biomass: {total_ERB:.2f} tCO2e")
        print(f"Average Extended Rotation Biomass per hectare: {average_ERB:.4f} tCO2e/ha")
    
    return total_results


## Step 5.3-ER: Define equation to calculate Total Benefits of ER
def calculate_total_benefits(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        benefit_live_biomass_ER = ERB

        total_benefits = benefit_live_biomass_ER
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "total_benefits": total_benefits
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Total benefits: {total_benefits:.2f} t CO2e")
    
    return total_results
  

### STEP 6-ER: Calculate total_benefits for bau and int and take difference between both scenarios
def calculate_all(scenario, log_level='info'):
    result = {}
    
    try:
        check_global_variables()
        
        steps = [
            ('average_agbd_tco2e', lambda: mean_annual_agc_stock(scenario=scenario, log_level=log_level)),
            ('mean_annual_tot_c_stock', lambda: mean_annual_tot_c_stock(result['average_agbd_tco2e'], scenario=scenario, log_level=log_level)),
            ('total_benefits', lambda: calculate_total_benefits(scenario=scenario, years=20, average_agbd_tco2e=result['average_agbd_tco2e'], log_level=log_level)),
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

####Calculate bau and intervention scenarios and difference between both

biomass = {
    "business_as_usual": calculate_all(scenario="business_as_usual"),
    "intervention": calculate_all(scenario="intervention")
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

    difference_ER = inter_item['annual_change'] - bau_item['annual_change']
    
    biomass_co2_result[aoi_id] = difference_ER
    biomass_co2_sd[aoi_id] = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96

print("Biomass CO2 Result ER:", biomass_co2_result)
print("Biomass CO2 Standard Deviation ER:", biomass_co2_sd)




########### STOP LOGGING (SL) SUB-INTERVENTION ###########
## Baseline: Conventional Logging Emissions
## Sub-Intervention: No emissions due to Stop Logging

## Step 4.1-SL: Define Baseline Conventional Logging Emissions equation
## Emissions (t CO2e) = Harvest area * (result from Step 4.1 + result from Step 4.2)
def calculate_conventional_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        forest_area = subregion["forest_area"]
        
        AHA_perc = float(scenario_data.get("AHA_perc", 0)) # this should be user input, as mentioned in User Input commented part of python code

        AHA = forest_area * AHA_perc ## edited this here to accomodate suggestion

    # Get results from Step 1-ER:
    subregion_results = mean_annual_agc_stock(scenario, log_level)
    total_c_stock_results = mean_annual_tot_c_stock(subregion_results, scenario, log_level)
    
    SL_emissions_results = []
    
    for c_stock in zip(total_c_stock_results):
        aoi_id = c_stock['aoi_id']
        forest_area = c_stock['area']
        
        # Get the harvest area from scenario data or use a default value
        AHA_perc = float(scenario_data.get("AHA_perc", area))
        
        # Calculate emissions
        total_c_stock = c_stock['total_carbon_stock']
        
        emissions = (AHA_perc * forest_area) * total_c_stock
        
        emissions_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "harvest_area": harvest_area,
            "total_carbon_stock": total_c_stock,
            "emissions": emissions
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Harvest Area: {harvest_area:.2f} ha")
            print(f"  Total Carbon Stock: {total_c_stock:.2f} tCO2e/ha")
            print(f"  Emissions from Conventional Logging: {emissions:.2f} tCO2e")
    
    return SL_emissions_results

## STEP 5-SL: FINAL CALCULATIONS TO OBTAIN CO2 IMPACT OF STOP LOGGING
def calculate_all(scenario, log_level='info'):
    result = {}
    
    try:
        check_global_variables()
        
        steps = [
            ('average_agbd_tco2e', lambda: mean_annual_agc_stock(scenario=scenario, log_level=log_level)),
            ('mean_annual_tot_c_stock', lambda: mean_annual_tot_c_stock(result['average_agbd_tco2e'], scenario=scenario, log_level=log_level)),
            ('SL_emissions_results', lambda: calculate_conventional_logging_emissions(scenario=scenario, years=20, average_agbd_tco2e=result['average_agbd_tco2e'], log_level=log_level)),
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

#### Step 6: Calculate bau and intervention scenarios and difference between both

biomass = {
    "business_as_usual": calculate_all(scenario="business_as_usual"),
    "intervention": calculate_all(scenario="intervention")
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

    difference_SL = inter_item['annual_change'] - bau_item['annual_change']
    
    biomass_co2_result[aoi_id] = difference_SL
    biomass_co2_sd[aoi_id] = (abs(difference_SL) * biomass_co2_result_error_positive / 100) / 1.96

print("Biomass CO2 Result SL:", biomass_co2_result)
print("Biomass CO2 Standard Deviation SL:", biomass_co2_sd)



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
    # NOTE: these are not required with the aoi_id proposed
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
