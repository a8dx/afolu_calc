# /***********************************************************************************************************
#   Copyright © 2024 Mathematica, Inc. This software was developed by Mathematica as part of the AFOLU GHG Calculator 
# project funded by USAID through Contract No. 51964. This code cannot be copied, distributed or used without 
# the express written permission of Mathematica, Inc.
# ***********************************************************************************************************/

# /**********************************************************************************************************
# Filename: Forest Conservation by Improved Cookstoves (FCIC)_test.py
# Author: Barbara Bomfim
# Date Started: 07/12/2024
# Last Edited: 07/15/2024
# Purpose: AFOLU GHG Calculations for Forest Conservation by Improved Cookstoves (FCIC) Interventions
# **********************************************************************************************************/

#### Required Inputs ####
# Climate zone and soil type: These are contained in the geography_mask_area_info json
#   file, which sources data from Dynamic World land use data,  data layers
# User inputs: User inputs depend on sub-intervention category
#   - Initial:
#           -land area: area under business-as-usual land use
#           -forest_type: forest type selection from dropdown list
#.          -forest_management_type: Natural, Plantation
#.  -Baseline data_FP:
#           -forest_area: Total forest area affected by fuelwood extraction (hectares)
#           -forest_type: forest type selection from dropdown list
#.          -forest_management_type: Natural, Plantation
#           -deforested_for_fuelwood: Percentage of harvested wood used for fuelwood
#           -fuelwood_extraction_rate: Annual rate of fuelwood extraction (hectares per year)
#           -R: fraction of rural users using cookstoves (Number of households using cookstoves relative to total households in rural area)
#           -fuel_type: wood, charcoal
#           -coosktove_performance: stove thermal efficiency (improved = 0.1)
#           -D: wood density (t/m3)
#           -fuel_consumption: Household fuel consumption (t/yr)
#           -kiln_yield_charcoal: Kiln yield for charcoal (1t of charcoal = Xt of wood)
#           -lfuelwood: Biomass removal rates for fuelwood (tC/ha)
#           -fNRB: Fraction of non-renewable biomass (Asia = 0.17, LAC = 0.33, Sub-Saharan Africa = 0.39, Global = 0.32)
#           -NCVbase = Net Caloric Value used in baseline
#           -nbase = thermal efficiency of baseline cookstove
#  -Intervention requires:
#           -forest_area: Total forest area affected by fuelwood extraction (hectares)
#           -forest_type: forest type selection from dropdown list
#.          -forest_management_type: Natural, Plantation
#           -R: fraction of rural users using cookstoves (Number of households using cookstoves relative to total households in rural area)
#           -number_improved_cookstoves: Number of households targeted with improved cookstoves
#           -displacement: percentage of household cooking displaced by improved stoves (percent)
#           -improved_cookstove_lifespan: Expected operation lifespan of improved cookstoves
#           -fuel_type: wood, charcoal
#           -coosktove_performance: stove thermal efficiency (improved = 0.35)
#           -D: wood density (t/m3)
#           -fuel_consumption: Household fuel consumption (t/yr)
#           -kiln_yield_charcoal: Kiln yield for charcoal (1t of charcoal = Xt of wood)
#           -lfuelwood: Biomass removal rates for fuelwood (tC/ha)
#           -fNRB: Fraction of non-renewable biomass (Asia = 0.17, LAC = 0.33, Sub-Saharan Africa = 0.39, Global = 0.32)
#           -NCVimpr = Net Caloric Value used in improved cookstove
#           -nimpr = thermal efficiency of improved cookstove

#### Parameters and Paths ####
# Calculations requires the following parameters and their associated uncertainty:

##Baselines and Improved Cookstoves
# "fNRBr"
# "fNRBu"
# "R":
# "Bbase"
# "NCVbase"
# "EFdirect_base"
# "EFkiln_base"
# "Ybase"

##N2O
# burning_n2o_ef: Emissions factor for N2O from burning (g N2O/kg dry matter burnt)
# combustion_factor: combustion factor for fire
# fuel_biomass: amount of biomass available for burning (t dry matter/ha)
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

###### Annual GHG impact GEDI data_FP and further calculations ######

## Import libraries
import json
import os

import pandas as pd
import numpy as np
from scipy.stats import norm
from concurrent.futures import ProcessPoolExecutor
# import ee
import re
import math

from helpersFCIC import convert_to_c, convert_to_co2e, convert_to_bgb, get_carbon_stock, sanitize_filename

## JSON data input #####

# JSON file path - Get the current working directory
current_directory = os.getcwd()
print("Current working directory:", current_directory)
# Define the relative file path
file_name = 'FCIC_test.json'
file_path = os.path.join(current_directory, file_name)
json_file_path = './FCIC_test.json'

# Load data_FP from the JSON file
with open(file_path, 'r') as file:
        data = json.load(file)

data 

### FCIC CO2 calculations ###

## conect from AGB data from GEE
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

#To have forest AOI carbon stock data - even though not needed in calculations below
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

####FCIC CALCULATIONS #####################
### Step 1: Calculate adjusted fNRB

#Define equation to calculate adjusted fNRB
def calculate_fNRBadj(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
    
        fNRBr = scenario_data.get("fNRBr")
        fNRBu = scenario_data.get("fNRBu")
        R = scenario_data.get("R")
    
        fNRBadj = fNRBr * R + fNRBu * (1 - R)
    
        total_results.append({
            "aoi_id": aoi_id,
            "area": area,
            "fNRBadj": fNRBadj
        })
        
        if log_level == 'debug':
            print(f"Adjusted fraction of fuelwood consumption for subregion {aoi_id}: {fNRBadj:.2f} %/yr")

    return total_results

#Define equation to calculate Non-renewable emissions (NRemissions)
## Baseline
def calculate_NRbase_emissions(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
    
        fNRBr = scenario_data.get("fNRBr")
        fNRBu = scenario_data.get("fNRBu")
        R = scenario_data.get("R")
        Bbase = scenario_data.get("Bbase")
        NCVbase = scenario_data.get("NCVbase")
        EFdirect_base = scenario_data.get("EFdirect_base")
        EFkiln_base = scenario_data.get("EFkiln_base")
        Ybase = scenario_data.get("Ybase")
    
        fNRBadj = fNRBr * R + fNRBu * (1 - R)
    
        edirect_base = Bbase * NCVbase * EFdirect_base * fNRBadj
    
        ekiln_base = Bbase * EFkiln_base * fNRBadj
    
        edamage_base = (Bbase / Ybase) * 0.32 * 0.47 * (44 / 12) * fNRBadj
    
        etotal_base = edirect_base + ekiln_base + edamage_base
    
        total_results.append({
        "fNRBadj": fNRBadj,
        "edirect_base": edirect_base,
        "ekiln_base": ekiln_base,
        "edamage_base": edamage_base,
        "etotal_base": etotal_base
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Total emissions baseline: {etotal_base:.2f} tCO2e/y")
    
    NRemissions_base = edirect_base + ekiln_base + edamage_base + etotal_base
    
    if log_level == 'debug':
        # print(f"Subregion {subregion['aoi_id']}:")
        print(f"Non-renewable cookstove emissions from one household: {NRemissions_base:.2f} tCO2e/y")
    
    return total_results


#Define equation to calculate the difference in thermal efficiency of baseline and improved stove (Bimpr)
def calculate_Bimpr(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
    
        Bbase = scenario_data.get("Bbase")
        NCVbase = scenario_data.get("NCVbase")
        NCVimpr = scenario_data.get("NCVimpr")
        nbase = scenario_data.get("nbase")
        nimpr = scenario_data.get("nimpr")
        
        Bimpr = Bbase * (nbase / nimpr) * (NCVbase / NCVimpr)
    
        total_results.append({
            "Bimpr": Bimpr
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Biomass consumption in improved cookstove: {Bimpr:.2f} t/y")
    
    return total_results

#Define equation to calculate Non-renewable emissions (NRemissions)
## Improved coockstoves
def calculate_NRimpr_emissions(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in total_results:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
 
        fNRBr = scenario_data.get("fNRBr")
        fNRBu = scenario_data.get("fNRBu")
        R = scenario_data.get("R")
        Bimpr = subregion["Bimpr"]
        NCVimpr = scenario_data.get("NCVimpr")
        EFdirect_impr = scenario_data.get("EFdirect_impr")
        EFkiln_impr = scenario_data.get("EFkiln_impr")
        Yimpr = scenario_data.get("Yimpr")
 
        fNRBadj = fNRBr * R + fNRBu * (1 - R)
    
        edirect_impr = Bimpr * NCVimpr * EFdirect_impr * fNRBadj
    
        ekiln_impr = Bimpr * EFkiln_impr * fNRBadj
    
        edamage_impr = (Bimpr / Yimpr) * 0.32 * 0.47 * (44 / 12) * fNRBadj
    
        etotal_impr = edirect_impr + ekiln_impr + edamage_impr
    
        total_results.append({
        "fNRBadj": fNRBadj,
        "edirect_impr": edirect_impr,
        "ekiln_impr": ekiln_impr,
        "edamage_impr": edamage_impr,
        "etotal_impr": etotal_impr
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Total emissions improved cookstove: {etotal_impr:.2f} tCO2e/y")
    
    NRemissions_impr = edirect_impr + ekiln_impr + edamage_impr + etotal_impr
    
    if log_level == 'debug':
        # print(f"Subregion {subregion['aoi_id']}:")
        print(f"Non-renewable cookstove emissions from one household: {NRemissions_impr:.2f} tCO2e/y")
    
    return total_results

#Define equation to calculate the annual savins in non-renewable CO2 emissions per household
def calculate_Bimpr(total_results, scenario, log_level='info'):
    global data
    
    scenario_data = data.get("scenarios", {}).get(scenario, [{}])[0]
    aoi_subregions = scenario_data.get("aoi_subregions", [])
    
    total_results = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
    
        Bbase = scenario_data.get("Bbase")
        NCVbase = scenario_data.get("NCVbase")
        NCVimpr = scenario_data.get("NCVimpr")
        nbase = scenario_data.get("nbase")
        nimpr = scenario_data.get("nimpr")
        
        esave = (etotal.base - etotal.impr) * Disp
    
        total_results.append({
            "Bimpr": Bimpr
        })
        
        if log_level == 'debug':
            print(f"Subregion {aoi_id}:")
            print(f"  Biomass consumption in improved cookstove: {Bimpr:.2f} t/y")
    
    return total_results





####### FOREST CONSERVATION BY IMPROVED COOKSTOVE INTERVENTION: Annual CO2 Emissions Calculations #######
#######FCIC Emissions ###############

def calculate_all_FCIC_bs(scenario):
    result_bs_FCIC = {}
    
    fNRBadj = calculate_fNRBadj(scenario=scenario)
    result_bs_FCIC["fNRBadj"] = fNRBadj
    
    etotal_base_d = calculate_baseline_emissions(scenario=scenario)
    result_bs_FCIC["etotal_base_d"] = etotal_base_d

    return result_bs_FCIC

bau_result_bs_FCIC = calculate_all_FCIC_bs(scenario="business_as_usual")
bau_result_bs_FCIC["etotal_base_d"]
print(f'Annual Baseline CO2e emissions (CO2_thayr): {bau_result_bs_FCIC["etotal_base_d"]}')


###### FCIC CO2 IMPACT  ########
### Calculate difference between intervention and business_as_usual

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

######## GHG CALCULATIONS #########

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
    result_df.to_csv('result.csv', index=False)
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




