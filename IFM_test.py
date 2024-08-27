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
#           -D: wood density (t/m3)
#           -VolExt_before = volume timber over bark extracted before project intervention (m3 ha-1)
#           -AHA (Annual Harvest Area): in ha - user input. 
#               -- Will attribute a value in the JSON but will indicate this in the code.
#           -C_timber: Proportion of total carbon extracted that resides in each wood product class 
#                 -- (Roundwood: Tropical "0.55", Temperate "0.45", Boreal "0.42", Sawnwood: Tropical 0.58, Temperate 0.48, Boreal 0.44, Woodbase panels: Tropical 0.62, Temperate 0.52, Boreal 0.52)
#           -prop_ox: Carbon emitted due to short-term oxidation of wood products "0.25" - the average across: Sawnwood 0.2; Woodbase panels 0.1; Other industrial roundwood 0.3; Paper and paperboard 0.4)
#           -prop_ox_long: Carbon in additional oxidized fraction "0.947" - the average across: Sawnwood 0.84; Woodbase panels 0.97; Other industrial roundwood 0.99; Paper and paperboard 0.99
#           -RL: Number of years in baseline rotation length
#           -GR: Growth rate
#           -BEF: Biomass Expansion Factor
#           -half_life_wp = Half-life of wood products
#  -Intervention requires:
#           -forest_area: forest area under intervention
#.          -forest_management_type: Natural, Plantation
#           -forest_type: forest type prior to intervention
#           -D: wood density (t/m3)
#           -VolExt_before = volume timber over bark extracted before project intervention (m3 ha-1)
#           -AHA (Annual Harvest Area): in ha - user input. 
#               -- Will attribute a value in the JSON but will indicate this in the code.
#           -C_timber: Proportion of total carbon extracted that resides in each wood product class 
#                 -- (Roundwood: Tropical "0.55", Temperate "0.45", Boreal "0.42", Sawnwood: Tropical 0.58, Temperate 0.48, Boreal 0.44, Woodbase panels: Tropical 0.62, Temperate 0.52, Boreal 0.52)
#           -prop_ox: Carbon emitted due to short-term oxidation of wood products "0.25" - the average across: Sawnwood 0.2; Woodbase panels 0.1; Other industrial roundwood 0.3; Paper and paperboard 0.4)
#           -prop_ox_long: Carbon in additional oxidized fraction "0.947" - the average across: Sawnwood 0.84; Woodbase panels 0.97; Other industrial roundwood 0.99; Paper and paperboard 0.99
#           -RL: Number of years in baseline rotation length
#           -GR: Growth rate
#           -BEF: Biomass Expansion Factor
#           -half_life_wp = Half-life of wood products

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

from helpersIFM import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock, sanitize_filename

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
file_name = 'IFM_test.json'
file_path = os.path.join(current_directory, file_name)
json_file_path = './IFM_test.json'
# Load data from the JSON file
with open(file_path, 'r') as file:
        data = json.load(file)

### STEP 2 - Calculate mean annual AGC stock in tCO2e/ha####

def mean_annual_agc_stock(scenario, log_level='info'):
    global data
    
    scenario_data = data["scenarios"][scenario][0]
    
    average_agbd_value = average_agbd_dict.get('agbd_mean', 0)
    area_converted = polygon.area().getInfo()/10000
    area_converted_yr = float(area_converted)
    
    CF = float(scenario_data.get("CF", 0))
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

### STEP 3 - Calculate mean annual Total C stock (sum AGC and BGC) ####

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

######## Annual CO2 Calculations for both Baseline and RIL intervention #########
## Equations:
# Logging emissions - baseline and RIL same equations
# Carbon stored in long-term wood products - baseline and RIL same equations
# Carbon emissions from incidental damage - baseline and RIL same equations

## Define Baseline Conventional Logging Emissions equation
def calculate_logging_emissions(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AHA = float(scenario_data.get("AHA", 0))
        VolExt = float(scenario_data.get("VolExt", 0))
        D = float(scenario_data.get("D", 0))
        ElE_factor_1 = float(scenario_data.get("ElE_factor_1", 0))
        ElE_factor_2 = float(scenario_data.get("ElE_factor_2", 0))
        
        ELE = (ElE_factor_1 * D) - ElE_factor_2
        TimberTree = (AHA * VolExt * ELE) * 44/12
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "TimberTree": TimberTree_conv,
            "ELE": ELE
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Extracted Log Emissions (ELE): {ELE:.4f} tC/m3 extracted")
            print(f"  Logging emissions : {TimberTree:.2f} tCO2e")
    
    return total_results

#Define equation to calculate carbon stored in long-term wood products
def calculate_carbon_wood_products(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        C_timber = float(scenario_data.get("C_timber", 0))
        r_w = float(scenario_data.get("r_w", 0))
        prop_ox = float(scenario_data.get("prop_ox", 0))
        prop_ox_long = float(scenario_data.get("prop_ox_long", 0))
        
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

#Define equation to calculate Carbon Emissions from Incidental Damage
def calculate_incidental_damage(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        TreeCarbonStock = float(scenario_data.get("TreeCarbonStock", 0))
        VolExt = float(scenario_data.get("VolExt", 0))
        LDF_factor_1 = float(scenario_data.get("LDF_factor_1", 0))
        LDF_factor_2 = float(scenario_data.get("LDF_factor_2", 0))
        #check how to get average_total_tco2e
        
        LDF = LDF_factor_1 * average_total_tco2e + LDF_factor_1
        IncidentalDamage = (LDF * VolExt) * 44/12
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "LDF": LDF,
            "IncidentalDamage": IncidentalDamage
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Tree Carbon Stock: {average_total_tco2e:.2f} tC/ha")
            print(f"  Volume Extracted Before: {VolExt:.2f} m3/ha")
            print(f"  Logging Damage Factor (LDF): {LDF:.4f}")
            print(f"  Incidental Damage Emissions: {IncidentalDamage:.2f} tCO2e")
    
    return total_results

#Define equation to calculate Conventional carbon emissions from Logging Infrastructure
def calculate_logging_infrastructure(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AHA = float(scenario_data.get("AHA", 0))
        VolExt_before = float(scenario_data.get("VolExt_before", 0))
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
            print(f"  Volume Timber Extracted (VolExt): {VolExt_before:.2f} m3/ha")
            print(f"  Skids emissions: {Skids:.2f} tC")
            print(f"  Roads and Decks Emissions: {RoadsDeck:.2f} tC")
            print(f"  Logging Infrastructure Emissions: {LoggingInfrastructure:.2f} tC")
    
    return total_results

# Estimate Total Logging Emissions
def calculate_total_logging_emissions(total_results, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    # Initialize the list to store the results
    total_results = []
    
    # Iterate through each result in total_results
    for result in total_results:
        aoi_id = result["aoi_id"]
        area = result["area"]
        
        # Retrieve precomputed values from the total_results
        TimberTree = float(result.get("TimberTree", 0))
        WoodProducts = float(result.get("WoodProducts", 0))
        IncidentalDamage = float(result.get("IncidentalDamage", 0))
        LoggingInfrastructure = float(result.get("LoggingInfrastructure", 0))
        
        # Calculate Emissions Conventional Logging
        Logging_Emissions = ((TimberTree - WoodProducts) + IncidentalDamage + LoggingInfrastructure) * (44 / 12)
        
        # Append the results
        emissions_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "TimberTree": TimberTree,
            "WoodProducts": WoodProducts,
            "IncidentalDamage": IncidentalDamage,
            "LoggingInfrastructure": LoggingInfrastructure,
            "Logging_Emissions": Logging_Emissions
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Timber Tree Emissions: {TimberTree:.2f} t C")
            print(f"  Carbon stored in wood products: {WoodProducts:.2f} t C")
            print(f"  Incidental Damage Emissions: {IncidentalDamage:.2f} t C")
            print(f"  Logging Infrastructure Emissions: {LoggingInfrastructure:.2f} t C")
            print(f"  Emissions Conventional Logging: {Logging_Emissions:.2f} tCO2e")
    
    return total_results

###### CALCULATE SUB-INTERVENTIONS CO2 IMPACT #######

#####Calculate all IFM - RIL variables ####
def calculate_all(scenario, log_level='info'):
    result = {}
    
    average_agbd_tco2e = mean_annual_agc_stock(scenario=scenario, log_level=log_level)
    result["average_agbd_tco2e"] = average_agbd_tco2e
    
    mean_annual_tot_c_stock_result = mean_annual_tot_c_stock(average_agbd_tco2e, scenario=scenario, log_level=log_level)
    result["mean_annual_tot_c_stock"] = mean_annual_tot_c_stock_result
    
    TimberTree = calculate_logging_emissions(scenario=scenario, log_level=log_level)
    result["TimberTree"] = TimberTree
    if log_level == 'debug':
        print(f"Annual Logging Emissions: {TimberTree} tCO2e/yr")
    
    WoodProducts = calculate_carbon_wood_products(scenario=scenario, log_level=log_level)
    result["WoodProducts"] = WoodProducts
    if log_level == 'debug':
        print(f"Carbon stored in wood products: {WoodProducts} tCO2e")
        
    IncidentalDamage = calculate_incidental_damage(scenario=scenario, log_level=log_level)
    result["IncidentalDamage"] = IncidentalDamage
    if log_level == 'debug':
        print(f"Logging Emissions due to Incidental Damage: {IncidentalDamage} tCO2e/yr")

    LoggingInfrastructure = calculate_logging_infrastructure(scenario=scenario, log_level=log_level)
    result["LoggingInfrastructure"] = LoggingInfrastructure
    if log_level == 'debug':
        print(f"Logging Emissions due to Incidental Damage: {LoggingInfrastructure} tCO2e/yr")

    Logging_Emissions = calculate_total_logging_emissions(scenario=scenario, log_level=log_level)
    result["Logging_Emissions"] = Logging_Emissions
    if log_level == 'debug':
        print(f"Logging Emissions due to Incidental Damage: {Logging_Emissions} tCO2e/yr")

    return result

biomass = {}
biomass["business_as_usual"] = calculate_all(scenario="business_as_usual")
biomass["intervention"] = calculate_all(scenario="intervention")

biomass_co2_result = {}

for inter_item, bau_item in zip(biomass["intervention"]["Logging_Emissions"], 
                                biomass["business_as_usual"]["Logging_Emissions"]):
    if inter_item['aoi_id'] != bau_item['aoi_id']:
        raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
    
    aoi_id = inter_item['aoi_id']
    difference = inter_item['Logging_Emissions'] - bau_item['Logging_Emissions']
    
    biomass_co2_result[aoi_id] = difference inter_result = calculate_all(scenario="intervention")

print(biomass_co2_result)

biomass_co2_result_error_positive = 10  # Assuming this is a percentage

biomass_co2_sd = {}

for aoi_id, difference in biomass_co2_result.items():
    # Calculate SD for each aoi_id
    sd = (abs(difference) * biomass_co2_result_error_positive / 100) / 1.96
    biomass_co2_sd[aoi_id] = sd

####### EXTENDED ROTATION INTERVENTION ##########

## Defining equations to be used if Extended Rotation is the selected sub-intervention

#Define equation to calculate baseline total tree carbon stock
def calculate_baseline_biomass_carbon_stock(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        average_agbd_tco2e = float(scenario_data.get("AGB", 0))
        ratio = float(scenario_data.get("ratio_below_ground_biomass_to_above_ground_biomass", 0))
        
        BL = average_agbd_tco2e * (1 + ratio)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "BL": BL
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Baseline live biomass: {BL:.2f} t CO2e")
    
    return total_results

#Define equation to calculate wood products - baseline

def calculate_wood_products_baseline(scenario, t=None, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AGB = float(scenario_data.get("AGB", 0))
        BEF = float(scenario_data.get("BEF", 0))
        ConvEff = float(scenario_data.get("ConvEff", 0))
        
        WPBL0 = AGB / BEF * ConvEff
        
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

#Define equation to calculate Biomass Accumulation due to ER
def calculate_extended_rotation_live_biomass(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        average_agbd_tco2e = float(scenario_data.get("AGB", 0))
        ERY = float(scenario_data.get("ERY", 0))
        GR = float(scenario_data.get("GR", 0))
        RS = float(scenario_data.get("RS", 0))
        
        ERB = (average_agbd_tco2e + (ERY * GR)) * (1 + RS)
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "ERB": ERB
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Extended rotation live biomass: {ERB:.2f} t CO2e")
    
    return total_results

#Define equation to calculate the the long term live tree biomass
def calculate_long_term_live_biomass(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        BL = float(scenario_data.get("BL", 0))
        ERB = float(scenario_data.get("ERB", 0))
        
        LTB = BL / 2
        LTERB = ERB / 2
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "LTB": LTB,
            "LTERB": LTERB
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Long-term biomass in baseline: {LTB:.2f} t CO2e")
            print(f"  Long-term biomass in extended rotation: {LTERB:.2f} t CO2e")
    
    return total_results

#Define equation to calculate carbon in wood products from extended rotation

import math

def calculate_wood_products_extended_rotation(scenario, t=None, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        AGB = float(scenario_data.get("AGB", 0))
        ERY = float(scenario_data.get("ERY", 0))
        GR = float(scenario_data.get("GR", 0))
        BEF = float(scenario_data.get("BEF", 0))
        
        WPER0 = (AGB + (ERY * GR)) / BEF
        
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
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        LTERB = float(scenario_data.get("LTERB", 0))
        LTB = float(scenario_data.get("LTB", 0))
        effectiveness = float(scenario_data.get("effectiveness", 0))
        
        benefit_live_biomass_ER = ((LTERB * area) - (LTB * area)) * (44 / 12) * effectiveness
        
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "benefit_live_biomass_ER": benefit_live_biomass_ER
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Benefit long term live biomass extended rotation: {benefit_live_biomass_ER:.2f} t CO2e")
    
    return total_results

def calculate_benefit_extended_rotation_wood_products(scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        LTWPER = float(scenario_data.get("LTWPER", 0))
        LTWPB = float(scenario_data.get("LTWPB", 0))
        effectiveness = float(scenario_data.get("effectiveness", 0))
        
        benefit_wood_products_ER = ((LTWPER * area) - (LTWPB * area)) * (44 / 12) * effectiveness
        
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

for inter_item, bau_item in zip(biomass_er["intervention"]["Logging_Emissions"], 
                                biomass_er["business_as_usual"]["Logging_Emissions"]):
    if inter_item['aoi_id'] != bau_item['aoi_id']:
        raise ValueError(f"Mismatched aoi_ids: {inter_item['aoi_id']} and {bau_item['aoi_id']}")
    
    aoi_id = inter_item['aoi_id']
    difference = inter_item['Logging_Emissions'] - bau_item['Logging_Emissions']
    
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
    # df['SOC_ref_tonnes_C_ha_sd'] = (df['high_activity_clay_soils_HAC_error_positive'] / 100 * df['SOC_ref_tonnes_C_ha']) / 1.96
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
        # # SOC Stock Change
        # SOCREF = np.random.normal(temp_bau['SOC_ref_tonnes_C_ha'].values[0], temp_bau['SOC_ref_tonnes_C_ha_sd'].values[0])
        # FLUbau = temp_bau['FLU'].values[0]
        # FMGbau = np.random.normal(temp_bau['FMG'].values[0], temp_bau['FMG_sd'].values[0])
        # FIbau = np.random.normal(temp_bau['FI'].values[0], temp_bau['FI_sd'].values[0])
        # FLUint = temp_int['FLU'].values[0]
        # FMGint = np.random.normal(temp_int['FMG'].values[0], temp_int['FMG_sd'].values[0])
        # FIint = np.random.normal(temp_int['FI'].values[0], temp_int['FI_sd'].values[0])
        # SOCbau = SOCREF * FLUbau * FMGbau * FIbau
        # SOCint = SOCREF * FLUint * FMGint * FIint
        # dSOC = SOCbau - SOCint
        # results_unc["SOC"].append(dSOC)
        # 
        # # here we add the biomass specific to the aoi_id
        # # NOTE: this number is too large check please
        # results_unc["totalC"].append(dSOC * 44/12 + biomass_co2_result[aoi_id]) 
        # 
        # # N2O emissions
        # # Calculate change in N sources
        # FSN = (temp_int['n_fertilizer_amount'].values[0] * temp_int['n_fertilizer_percent'].values[0] / 100 -
        #        temp_bau['n_fertilizer_amount'].values[0] * temp_bau['n_fertilizer_percent'].values[0] / 100) if "nutrient management" in intervention_subcategory else 0
        # FON = (temp_int['org_amend_rate'].values[0] * temp_int['org_amend_npercent'].values[0] / 100 -
        #        temp_bau['org_amend_rate'].values[0] * temp_bau['org_amend_npercent'].values[0] / 100) if "nutrient management" in intervention_subcategory else 0
        # FSOM = dSOC / 10 * 1000  # FSOM
        # if "change in livestock type or stocking rate" in intervention_subcategory:
        #     FPRP = (temp_int['live_weight'].values[0] * temp_int['Nex'].values[0] * temp_int['stocking_rate'].values[0] -
        #             temp_bau['live_weight'].values[0] * temp_bau['Nex'].values[0] * temp_bau['stocking_rate'].values[0]) / 1000 * 365
        # else:
        #     FPRP = 0
        # 
        # # Calculate N emissions
        # EFdir = np.random.normal(temp_bau['ef_direct_n2o'].values[0], temp_bau['ef_direct_n2o_sd'].values[0])
        # FRAC_GASF = np.random.normal(temp_bau['frac_syn_fert'].values[0], temp_bau['frac_syn_fert_sd'].values[0])
        # FRAC_GASM = np.random.normal(temp_bau['frac_org_fert'].values[0], temp_bau['frac_org_fert_sd'].values[0])
        # EF_vol = np.random.normal(temp_bau['ef_vol'].values[0], temp_bau['ef_vol_sd'].values[0])
        # FRAC_LEACH = np.random.normal(temp_bau['frac_leach'].values[0], temp_bau['frac_leach_sd'].values[0])
        # EF_leach = np.random.normal(temp_bau['ef_leaching'].values[0], temp_bau['ef_leaching_sd'].values[0])
        # EF_PRP = np.random.normal(temp_bau['ef_prp'].values[0], temp_bau['ef_prp_sd'].values[0])
        # dirN = (FSN + FON + FSOM) * EFdir + FPRP * EF_PRP  # Direct emissions
        # volN = (FSN * FRAC_GASF + (FON + FPRP) * FRAC_GASM) * EF_vol  # Indirect volatilization
        # leachN = (FSN + FON + FSOM + FPRP) * FRAC_LEACH * EF_leach  # indirect leaching
        # N2O = (dirN + volN + leachN) / 1000 * 44 / 28  # sum and convert to tN2O
        # 
        # results_unc[f"N2O_{aoi_id}"].append(N2O)
        # 
        #FIXME: need to incorporeate fire management 
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

      

    # column_names = ['AOI', 'Area', 'Rep', 'SOC', 'totalC', f'N2O_y{aoi_id}', f'CH4_y{aoi_id}']
        results_unc_df = pd.DataFrame(results_unc)
    return results_unc_df


def generate_output():
    global biomass_co2_result
    json_data, df, AOIs, intervention_subcategory = load_json_data('IFM_test.json')
    
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
    result_df.to_csv('output\IFM_test_result.csv', index=False)
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
