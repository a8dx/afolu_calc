# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: sample_ARR_gee-agb_GHG-calculations_updated.py
# Author: Barbara Bomfim
# Date Started: 06/28/2024
# Last Edited: 07/11/2024
# Purpose: AFOLU GHG Calculations for Afforestation/Reforestation (A/R) Interventions
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data from Dynamic World land use data,  data layers
# User inputs: User inputs depend on sub-intervention category (planted or natural forest; 
#  tillage; nutrient management; fire management). 
#   -Baseline data:
#           -land area: area under business-as-usual land use
#           -land_use_type: previous land use type (Agroforestry, Annual cropland, Annual fallow, Degraded land, Flooded rice, Grassland)
#           -forest_type: selection from dropdown list
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#           -fire_used: yes, no
#  -All sub-interventions require:
#           -forest_area: area that is being changed from previous land use to new A/R land use
#           -number of improvements over initial scenario (0, 1 or 2)
#.    -For Forest creation sub-intervention:
#           -forest_creation: planted or natural forest (native forest or mangrove forest)
#           -forest_type: selection from dropdown list
#     -For tillage sub-intervention:
#           -tillage_type: Full Till, Reduced Till, No Till, Unknown
#.    -For nutrient management sub-intervention:
#           -ag_inputs: Low, Medium, High without manure, High with manure, Unknown
#     -For fire management, user must provide:
#.          -fire_used: yes, no
#           -frequency of burn (number of years)

#### Parameters and Paths ####
# Calculations requires the following parameters and their associated uncertainty:

##AGB
#AGBREF:AGB for the climate zone and soil type (tC/ha)

##BGB
#BGBREF:BGB for the climate zone and soil type (tC/ha)

##SOC
# SOCREF: Reference soil stock for the climate zone and soil type (t C/ha)
# FLU: Land use emissions factor (1 for planted forest, 1 for natural forest)
# FMG: Management factor for both business as usual and intervention scenario
# FMG: Tillage factor for both business as usual and intervention scenario (if 
  #    tillage is not a subintervention then FMG = 1)
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

## Annual GHG impact GEDI data and further calculations ###

## Import libraries
import json
import os

import pandas as pd
import numpy as np
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor
import ee
import re

from helpers import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock

### STEP 1 - Estimate mean annual AGB (in Mg/ha) using GEDI dataset ####
# GEE path: https://code.earthengine.google.com/?scriptPath=users%2Fbabomfimf%2FAGB-GEDI-L4A%3Atest-3_AGB_annual_mean
## Anthony's GEE: https://code.earthengine.google.com/44d6aaa5db4764b5b5f3825baf900d04

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
file_name = 'ARR_final_updated.json'
file_path = os.path.join(current_directory, file_name)
json_file_path = './ARR_final_updated.json'
# Load data from the JSON file
with open(file_path, 'r') as file:
        data = json.load(file)

### STEP 2 - Calculate mean annual AGC stock in tCO2e/ha####

def mean_annual_agc_stock(scenario):
    global data
    # Convert biomass stock to carbon stock in tCO2e/ha
    average_agbd_value = average_agbd_dict.get('agbd_mean',0)
    area_converted = print(f"{polygon.area().getInfo()/10000} ha")
    area_converted_dict = polygon.area().getInfo()/10000
    area_converted_yr = float(area_converted_dict)

    print(f"Average Aboveground Biomass Density: {average_agbd_value} Mg/ha")
    CF = data.get("scenarios").get(scenario)[0]["CF"]
    average_agbd_cstock = convert_to_c(average_agbd_value, CF)
    average_agbd_tco2e = convert_to_co2e(average_agbd_cstock)
    print(f"Average Aboveground Carbon Stock: {average_agbd_tco2e} tCO2e/ha")
    return average_agbd_tco2e


### STEP 3 - Calculate mean annual Total C stock (sum AGC and BGC) ####

def mean_annual_tot_c_stock(average_agbd_tco2e, scenario):
    global data
    # Calculate belowground carbon stock
    ratio = data.get("scenarios").get(scenario)[0]["ratio_below_ground_biomass_to_above_ground_biomass"]
    average_bgbd_tco2e = convert_to_bgb(average_agbd_tco2e, ratio)
    average_total_tco2e = average_agbd_tco2e + average_bgbd_tco2e# in tCO2e/ha

### STEP 4 - Calculate Annual CO2 impact (∆CG)####

#Annual-CO2-impact = [[Carbon-i - Carbon-0] - [Carbon-bau - DeltaCarbon-0]]/years
#Carbon-0: initial carbon stock (tC/ha) - average_agbd_cstock from above for year 0
#Carbon-i: carbon stock after intervention (tC/ha) - average_agbd_cstock for year i after intervention
#Carbon-bau: carbon stock under business-as-usual (without intervention) - average_agbd_cstock for year i without intervention
# Define equation to calculate Annual CO2 impact (∆CG)
def calculate_annual_co2_impact(scenario, years):
    """
    Calculate the annual CO2 impact using the provided equation.

    Parameters:
    scenario (str): Scenario in question
    years (int): Number of years over which the change is calculated.

    Returns:
    float: Annual CO2 impact in tCO2e/yr.
    """
    global data
    # Carbon stock at year i in tC/ha.
    carbon_i = get_carbon_stock(data, scenario=scenario, carbon_time_id="carbon_i")
    # Initial carbon stock in tC/ha.
    carbon_0 = get_carbon_stock(data, scenario=scenario, carbon_time_id="carbon_0")
    # Business-as-usual carbon stock in tC/ha.
    carbon_bau = get_carbon_stock(data, scenario=scenario, carbon_time_id="carbon_bau")
    # Change in carbon stock in the baseline scenario in tC/ha.
    delta_carbon_0 = data.get("scenarios").get(scenario)[0]["delta_carbon_0"]
    
    carbon_i_co2e = convert_to_co2e(carbon_i)
    carbon_0_co2e = convert_to_co2e(carbon_0)
    carbon_bau_co2e = convert_to_co2e(carbon_bau)
    delta_carbon_0_co2e = convert_to_co2e(delta_carbon_0)

    delta_CG = ((carbon_i_co2e - carbon_0_co2e) - (carbon_bau_co2e - delta_carbon_0_co2e)) / years
    print(f"Annual CO2 impact (∆CG): {delta_CG} tCO2e/yr in {years} years of intervention")
    return delta_CG

### Step 5 #### Estimate Initial change in biomass carbon stocks on converted land
#Adapt Equation 2.16 to estimate ∆C-CONVERSION as the immediate change in biomass due to conversion.
#∆C-CONVERSION = (AGB month after - AGB before) x area of land converted
def estimate_co2_conversion(scenario):
    """
    Estimate ∆C-CONVERSION as the immediate change in biomass due to conversion.

    Parameters:
    scenario (str): Scenario in question
    
    Returns:
    float: Immediate change in total biomass carbon due to conversion (tCO2e/yr).
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    # Aboveground biomass up to one month after the conversion (tC/ha).
    agb_month_after = scenario_data.get("agb_month_after")
    agb_month_after_co2e = convert_to_co2e(agb_month_after)

    # Aboveground biomass before conversion (tC/ha).
    agb_before = scenario_data.get("agb_before")
    agb_before_co2e = convert_to_co2e(agb_before)

    # ratio_below_ground_biomass_to_above_ground_biomass (dimensionless).
    ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")

    # Area of land converted in certain year (ha/yr).
    area_converted_yr = scenario_data.get("area_converted_yr")

    delta_co2_conversion = ((agb_month_after_co2e + (agb_month_after_co2e * ratio)) - (agb_before_co2e + (agb_before_co2e * ratio))) * area_converted_yr
    return delta_co2_conversion

### Step 6 #### Calculate Annual Decrease in Biomass Stock ####
#Esimate ∆CL using Equation 2.11: ∆CL = Lwood −removals + Lfuelwood + Ldisturbance
# Lwood-removals = annual carbon loss due to wood removals, tonnes C yr-1 (Equation 2.12) 
# Lfuelwood = annual biomass carbon loss due to fuelwood removals, tC/yr (Equation 2.13)
# Ldisturbance = annual biomass carbon losses due to disturbances, tonnes C yr-1 (See Equation 2.14)

#Defining Lwood-removals equation
def calculate_lwood_removals(scenario):
    """
    Calculate Lwood-removals using the provided equation.

    Parameters:
    scenario (str): Scenario in question

    Returns:
    float: Annual carbon loss due to wood removals in tCO2e/yr.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]

    # Annual wood removal in m3/yr.
    H = scenario_data.get("H")
    # Biomass Conversion and Expansion Factor from Table 4.5 in t d.m. removal/m3 removals.
    BCEFr = scenario_data.get("BCEFr")
    # Root-to-shoot ratio.
    ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
    # Carbon fraction.
    CF = scenario_data.get("CF")

    if None in [H, BCEFr, ratio, CF]:
        print(f"Missing data for scenario {scenario}. Please check your data.")
        return None


    return H * BCEFr * (1 + ratio) * CF * (44 / 12)

# Lfuelwood equation
def calculate_lfuelwood(scenario):
    """
    Calculate Lfuelwood using the provided equation.

    Parameters:
    scenario (str): Scenario in question

    Returns:
    float: Annual biomass carbon loss due to fuelwood removals in tCO2e/yr.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    
    # Annual volume of fuelwood removal of whole trees in m3/yr.
    FGtrees = scenario_data.get("FGtrees")
    # Annual volume of fuelwood removal of tree parts in m3/yr.
    FGpart = scenario_data.get("FGpart")
    # Biomass Conversion and Expansion Factor from Table 4.5 in t biomass removal/m3 removals.
    BCEFr = scenario_data.get("BCEFr")
    # Root-to-shoot ratio.
    ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
    # Wood density.
    D = scenario_data.get("D")
    # Carbon fraction.
    CF = scenario_data.get("CF")

    # Handling missing data
    if None in [FGtrees, FGpart, BCEFr, ratio, D, CF]:
        print(f"Missing data for calculating Lfuelwood. Please check your data for scenario {scenario}.")
        return None

    # Calculating Lfuelwood based on the provided equation
    Lfuelwood = ((FGtrees * BCEFr * (1 + ratio)) + FGpart * D) * CF * (44 / 12)

    return Lfuelwood

# Ldisturbance equation
def calculate_ldisturbance(scenario):
    """
    Calculate Ldisturbance using the provided equation.

    Parameters:
    scenario (str): Scenario in question

    Returns:
    float: Annual biomass carbon losses due to disturbances in tCO2e/yr.
    """
    global data
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]

    # Area affected by disturbances in ha/yr.
    Adisturbance = scenario_data.get("Adisturbance")
    
    # Average AGB of land areas affected by disturbance in tC/ha.
    Bw = scenario_data.get("Bw")
    
    # Root-to-shoot ratio (dimensionless).
    ratio = scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass")
    
    # Fraction of biomass carbon lost in disturbance (dimensionless).
    fd = scenario_data.get("fd")

    # Handling missing data
    if None in [Adisturbance, Bw, ratio, fd]:
        print("Missing data for calculating Ldisturbance. Please check your data.")
        return None

    # Calculating Ldisturbance based on the provided equation
    Ldisturbance = Adisturbance * Bw * (1 + ratio) * fd * (44 / 12)

    return Ldisturbance

#∆CL equation
def calculate_delta_cl(lwood_removals, lfuelwood, ldisturbance):
    """
    Calculate the change in biomass stock (∆CL).

    Parameters:
    scenario (str): Scenario in question
    lwood_removals (float): Annual carbon loss due to wood removals in tCO2e/yr.
    lfuelwood (float): Annual biomass carbon loss due to fuelwood removals in tCO2e/yr.
    ldisturbance (float): Annual biomass carbon losses due to disturbances in tCO2e/yr.


    Returns:
    float: Change in biomass stock (∆CL) in tCO2e/yr.
    """
    # Calculating ΔCL based on the provided equation
    return lwood_removals + lfuelwood + ldisturbance  


### Step 7 #### Calculate annual change in carbon stocks in biomass on land converted to other land-use category
# Estimate ∆CB using Equation 2.15: ∆CB = ∆CG + ∆C-CONVERSION − ∆CL
#∆CB = annual change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr)
#∆CG = annual increase in carbon stocks in biomass due to growth on land converted to another land-use category (tCO2e/yr)
#∆C-CONVERSION = initial change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr)
#∆CL = annual decrease in biomass carbon stocks due to losses from harvesting, fuel wood gathering and disturbances on land converted to other land-use category (tCO2e/yr)

#∆CB Equation
def calculate_annual_change_in_carbon_stocks(delta_CG, delta_co2_conversion, delta_cl):
    """
    Calculate the annual change in carbon stocks in biomass on land converted to other land-use category.

    Parameters:
    delta_CG (float): Annual increase in carbon stocks in biomass due to growth on land converted to another land-use category (tCO2e/yr)
    delta_co2_conversion (float): Initial change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr)
    delta_cl (float): Annual decrease in biomass carbon stocks due to losses from harvesting, fuel wood gathering and disturbances on land converted to other land-use category (tCO2e/yr)

    Returns:
    float: Annual change in carbon stocks in biomass on land converted to other land-use category (tCO2e/yr)
    """
    return delta_CG + delta_co2_conversion - delta_cl


def calculate_all(scenario):
    result = {}
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario)
    result["average_agbd_tco2e"] = average_agbd_tco2e

    mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock

    delta_CG = calculate_annual_co2_impact(scenario=scenario, years=30)
    result["delta_CG"] = delta_CG

    delta_co2_conversion = estimate_co2_conversion(scenario=scenario)
    result["delta_co2_conversion"] = delta_co2_conversion

    # Calculate Lwood-removals
    lwood_removals = calculate_lwood_removals(scenario=scenario)
    result["lwood_removals"] = lwood_removals
    print(f"Annual biomass carbon loss due to wood removals (Lwood-removals): {lwood_removals} tCO2e/yr")

    # Calculate Lfuelwood
    lfuelwood = calculate_lfuelwood(scenario=scenario)
    result["lfuelwood"] = lfuelwood
    print(f"Annual biomass carbon loss due to fuelwood removals (Lfuelwood): {lfuelwood} tCO2e/yr")

    # Calculate Ldisturbance
    ldisturbance = calculate_ldisturbance(scenario=scenario)
    result["ldisturbance"] = ldisturbance
    print(f"Annual biomass carbon losses due to disturbances (Ldisturbance): {ldisturbance} tCO2e/yr")

    # Calculate ∆CL
    delta_cl = calculate_delta_cl(lwood_removals, lfuelwood, ldisturbance)
    result["delta_cl"] = delta_cl
    print(f"Biomass stock Loss (∆CL): {delta_cl} tCO2e/yr")

    # Calculate the annual change in carbon stocks
    annual_change_in_carbon_stocks = calculate_annual_change_in_carbon_stocks(delta_CG, delta_co2_conversion, delta_cl)
    result["annual_change_in_carbon_stocks"] = annual_change_in_carbon_stocks
    print(f"Annual change in total biomass carbon stock (∆CB): {annual_change_in_carbon_stocks} tCO2e/yr")
    
    return result

bau_result = calculate_all(scenario="business_as_usual")
print("\n\n")
inter_result = calculate_all(scenario="intervention")

#### Final biomass CO2 impact ####
bau_bio_co2 = bau_result["annual_change_in_carbon_stocks"]
inter_bio_co2 = inter_result["annual_change_in_carbon_stocks"]
ghg_result = inter_bio_co2 - bau_bio_co2
ghg_result

# Load JSON data from a file
with open('ARR_final_updated.json', 'r') as f:
    data = json.load(f)

# Extract the "years" variable from the "business_as_usual" and "intervention" scenarios
years_bau = data["scenarios"]["business_as_usual"][0]["years"]
years_intervention = data["scenarios"]["intervention"][0]["years"]

print("Years (Business As Usual):", years_bau)
print("Years (Intervention):", years_intervention)

print(f"The years value is: {years}")
print(f"Annual biomass CO2 impact (∆CB): {ghg_result} tCO2e/yr in {years_bau} years of intervention")

##### SOIL GHG CALCULATIONS ######

## Getting json data input ##
# Load data from JSON input
with open('ARR_final_updated.json') as f:
    data = json.load(f)

# Extract intervention parameters
interv_sub = data['intervention_subcategory']  # intervention type

# Convert common scenarios to DataFrame
common = pd.json_normalize(data['scenarios']['common'])

# Convert business-as-usual scenario to DataFrame
bau = pd.json_normalize(data['scenarios']['business_as_usual'])
bau = bau.drop(columns=['aoi_subregions'])
bau

# Extract AOI subregions for business-as-usual
bau_aoi = pd.json_normalize(data['scenarios']['business_as_usual'], 'aoi_subregions')
bau_aoi['aoi_id'] = range(1, len(bau_aoi) + 1)
bau_aoi

# Combine BAU and AOI subregions
bau = pd.concat([bau, bau], ignore_index=True)
bau = pd.concat([bau, bau_aoi], axis=1)
bau['scenario'] = 'business-as-usual'
bau

# Convert intervention scenario to DataFrame
int_df = pd.json_normalize(data['scenarios']['intervention'])
int_df = int_df.drop(columns=['aoi_subregions'])
int_df

# Extract AOI subregions for intervention
int_aoi = pd.json_normalize(data['scenarios']['intervention'], 'aoi_subregions')
int_aoi['aoi_id'] = range(1, len(int_aoi) + 1)
int_aoi

# Combine intervention and AOI subregions
int_df = pd.concat([int_df, int_df], ignore_index=True)
int_df = pd.concat([int_df, int_aoi], axis=1)
int_df['scenario'] = 'intervention'
int_df

# Combine BAU and intervention DataFrames
df = pd.concat([bau, int_df], ignore_index=True)

# Repeat common scenarios DataFrame for each entry
common_repeated = pd.concat([common] * len(df), ignore_index=True)
df = pd.concat([common_repeated, df], axis=1)

# Convert data types appropriately
df = df.apply(pd.to_numeric, errors='ignore')

# Extract unique AOI identifiers
AOIs = df['aoi_id'].unique()
AOIs

# Convert uncertainty values to standard deviation (sd)
# Uncertainty presented as 95% CI
unc_low = [col for col in df.columns if re.search("uncertainty_low", col)]
unc_low_name = df[unc_low].columns
unc_cols = [re.sub("uncertainty_low.*", "", col) for col in unc_low_name]

for a in unc_cols:
    df[f"{a}sd"] = (df[f"{a}uncertainty_high"] - df[f"{a}uncertainty_low"]) / 3.92

# Drop columns with "uncertainty" in their name
df = df[df.columns.drop(list(df.filter(regex='uncertainty')))]

# Uncertainty presented as 95% CI as percentage of mean
df["SOC_ref_tonnes_C_ha_sd"] = (df["low_activity_clay_soils_LAC_error_positive"] / 100 * df["SOC_ref_tonnes_C_ha"]) / 1.96

# Uncertainty presented as 2 sd as percentage of mean
df["FMG_sd"] = (df["FMG_error_positive"] / 100 * df["FMG"]) / 2
df["FI_sd"] = (df["FI_error_positive"] / 100 * df["FI"]) / 2

# Remove other uncertainty columns
df = df[df.columns.drop(list(df.filter(regex='error')))]

# Display structure of the DataFrame
print(df.info())

## Soil GHG Calculations ##
nx=500 # number of Monte Carlo iterations

def GHGcalc(i, df, nx):
    temp_bau = df[(df['aoi_id'] == i) & (df['scenario'] == 'business-as-usual')]
    temp_int = df[(df['aoi_id'] == i) & (df['scenario'] == 'intervention')]
    results_unc = np.zeros((nx, 44))
    
    results_unc[:, 0] = i
    results_unc[:, 1] = temp_bau['area'].values[0]
    results_unc[:, 2] = np.arange(1, nx + 1)

    for m in range(nx):
        ####SOC Stock Change####
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
        results_unc[m, 3] = (dSOC * 44 / 12) + ghg_result

        ####N2O emissions####
        # Calculate change in N sources
        FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
               temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100) if 'nutrient management' in interv_sub else 0
        FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
               temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100) if 'nutrient management' in interv_sub else 0
        FSOM = dSOC / 10 * 1000

        # Calculate N emissions
        EFdir = np.random.normal(temp_bau['ef_direct_n2o'].values[0], temp_bau['ef_direct_n2o_sd'].values[0])
        FRAC_GASF = np.random.normal(temp_bau['frac_syn_fert'].values[0], temp_bau['frac_syn_fert_sd'].values[0])
        FRAC_GASM = np.random.normal(temp_bau['frac_org_fert'].values[0], temp_bau['frac_org_fert_sd'].values[0])
        EF_vol = np.random.normal(temp_bau['ef_vol'].values[0], temp_bau['ef_vol_sd'].values[0])
        FRAC_LEACH = np.random.normal(temp_bau['frac_leach'].values[0], temp_bau['frac_leach_sd'].values[0])
        EF_leach = np.random.normal(temp_bau['ef_leaching'].values[0], temp_bau['ef_leaching_sd'].values[0])
        #EF_PRP = np.random.normal(temp_bau['ef_prp'].values[0], temp_bau['ef_prp_sd'].values[0])
        dirN = (FSN + FON + FSOM) * EFdir
        volN = (FSN * FRAC_GASF + (FON) * FRAC_GASM) * EF_vol
        leachN = (FSN + FON + FSOM) * FRAC_LEACH * EF_leach
        N2O = (dirN + volN + leachN) / 1000 * 44 / 28
        results_unc[m, 4:24] = N2O

        # Add change in N2O due to fire management
        if 'change in fire management' in interv_sub:
            fire_n2o_ef = np.random.normal(temp_bau['burning_n2o_ef_mean'].values[0], temp_bau['burning_n2o_ef_sd'].values[0])
            if temp_bau['fire_used'].values[0] == 'True':
                CF_bau = np.random.normal(temp_bau['combustion_factor_mean'].values[0], temp_bau['combustion_factor_sd'].values[0])
                MB_bau = np.random.normal(temp_bau['fuel_biomass_mean'].values[0], temp_bau['fuel_biomass_sd'].values[0])
                fireN2O_bau = MB_bau * CF_bau * fire_n2o_ef / 1000
                results_unc[m, 4] = N2O - fireN2O_bau
                fire_per_bau = temp_bau['fire_management_years'].values[0]
                fire_yrs_bau = []
                for y in range(1, 20):
                    if y % fire_per_bau == 0:
                        fire_yrs_bau.append(y)
                for y in fire_yrs_bau:
                    results_unc[m, 4 + y] = N2O - fireN2O_bau

            if temp_int['fire_used'].values[0] == 'True':
                CF_int = np.random.normal(temp_int['combustion_factor_mean'].values[0], temp_int['combustion_factor_sd'].values[0])
                MB_int = np.random.normal(temp_int['fuel_biomass_mean'].values[0], temp_int['fuel_biomass_sd'].values[0])
                fireN2O_int = MB_int * CF_int * fire_n2o_ef / 1000
                results_unc[m, 4] += fireN2O_int
                fire_per_int = temp_int['fire_management_years'].values[0]
                fire_yrs_int = []
                for y in range(1, 20):
                    if y % fire_per_int == 0:
                        fire_yrs_int.append(y)
                for y in fire_yrs_int:
                    results_unc[m, 4 + y] += fireN2O_int

        ####CH4 emissions####

        # Emissions from fire management
        if 'change in fire management' in interv_sub:
            fire_ch4_ef = np.random.normal(temp_bau['burning_ch4_ef_mean'].values[0], temp_bau['burning_ch4_ef_sd'].values[0])
            if temp_bau['fire_used'].values[0] == 'True':
                fireCH4_bau = MB_bau * CF_bau * fire_ch4_ef / 1000
                results_unc[m, 24] = CH4_live - fireCH4_bau
                for y in fire_yrs_bau:
                    results_unc[m, 24 + y] = CH4_live - fireCH4_bau

            if temp_int['fire_used'].values[0] == 'True':
                fireCH4_int = MB_int * CF_int * fire_ch4_ef / 1000
                results_unc[m, 24] += fireCH4_int
                for y in fire_yrs_int:
                    results_unc[m, 24 + y] += fireCH4_int

    results_unc_df = pd.DataFrame(results_unc, columns=['AOI', 'Area', 'Rep', 'SOC'] + [f'N2O_y{y}' for y in range(1, 21)] + [f'CH4_y{y}' for y in range(1, 21)])
    results_unc_df = results_unc_df.apply(pd.to_numeric, errors='coerce')
    return results_unc_df

# Extract mean and standard deviation from each subAOI and store in results dataframe
for i in range(len(AOIs)):
    mc.append(pd.DataFrame(np.random.randn(20, 44), columns=['Area'] + ['V'+str(i) for i in range(1, 44)]))

# Extract mean and sd from each subAOI and store in results df
resdf = np.zeros((len(AOIs), 166))
resdf[:, 0] = AOIs
columns = ['AOI', 'Area', 'CO2_thayr'] + \
          [f'N2O_thayr_y{i}' for i in range(1, 21)] + \
          [f'CH4_thayr_y{i}' for i in range(1, 21)] + \
          ['sdCO2ha'] + [f'sdN2Oha_y{i}' for i in range(1, 21)] + \
          [f'sdCH4ha_y{i}' for i in range(1, 21)] + \
          ['CO2_tyr'] + [f'N2O_tyr_y{i}' for i in range(1, 21)] + \
          [f'CH4_tyr_y{i}' for i in range(1, 21)] + \
          ['sdCO2'] + [f'sdN2O_y{i}' for i in range(1, 21)] + \
          [f'sdCH4_y{i}' for i in range(1, 21)]
resdf = pd.DataFrame(resdf, columns=columns)
resdf

for a in range(len(AOIs)):
    temp_res = mc[a]
    resdf.at[a, 'Area'] = temp_res['Area'].mean()
    for z in range(4, 44):
        resdf.at[a, columns[z-1]] = temp_res.iloc[:, z].mean()
        resdf.at[a, columns[z+40]] = temp_res.iloc[:, z].std()

for k in range(3, 84):
    resdf.iloc[:, 82+k] = resdf['Area'] * resdf.iloc[:, k]

resdf = resdf.drop(columns=['Area'])

# Convert to JSON output
restib = []

for index, row in resdf.iterrows():
    aoi_tib = {
        "AOI": f"Sub_AOI-{index+1}",
        "CO2_thayr": [row['CO2_thayr']] * 20,
        "sdCO2ha": [row['sdCO2ha']] * 20,
        "CO2_tyr": [row['CO2_tyr']] * 20,
        "sdCO2": [row['sdCO2']] * 20,
        "N2O_thayr": [row[f'N2O_thayr_y{y}'] for y in range(1, 21)],
        "sdN2Oha": [row[f'sdN2Oha_y{y}'] for y in range(1, 21)],
        "N2O_tyr": [row[f'N2O_tyr_y{y}'] for y in range(1, 21)],
        "sdN2O": [row[f'sdN2O_y{y}'] for y in range(1, 21)],
        "CH4_thayr": [row[f'CH4_thayr_y{y}'] for y in range(1, 21)],
        "sdCH4ha": [row[f'sdCH4ha_y{y}'] for y in range(1, 21)],
        "CH4_tyr": [row[f'CH4_tyr_y{y}'] for y in range(1, 21)],
        "sdCH4": [row[f'sdCH4_y{y}'] for y in range(1, 21)]
    }
    restib.append(aoi_tib)

json_output = json.dumps(restib, indent=4)
with open(f"/Users/barbarabomfim/Dropbox/R/afolu_calc/output/{interv_sub}.json", "w") as f:
    f.write(json_output)
    
json_output
