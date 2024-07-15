import re
# Equation to convert carbon stock in tC into tCO2e
def convert_to_c(biomass_stock, CF):
    """
    Convert biomass stock in Mg/ha to tC/ha using carbon fraction (CF).

    Parameters:
    biomass_stock (float): Biomass stock in Mg/ha.
    CF: Carbon Fraction

    Returns:
    float: Carbon stock in tC/ha.
    """
    return biomass_stock * CF

# Define equation to convert tCha to tCO2e/ha 
def convert_to_co2e(carbon_stock):
    """
    Convert carbon stock in tC/ha to CO2 equivalents in tCO2e/ha.

    Parameters:
    carbon_stock (float): Carbon stock in tC/ha.

    Returns:
    float: CO2 equivalents in tCO2e/ha.
    """
    return carbon_stock * (44 / 12)

# Define equation to convert tCha to tCO2e/ha 
def convert_to_bgb(abg_b_or_c, ratio_below_ground_biomass_to_above_ground_biomass):
    """
    Convert carbon stock in tC/ha to CO2 equivalents in tCO2e/ha.

    Parameters:
    abg_b_or_c (float): Aboveground biomass (Mg/ha) or Aboveground carbon stock in tC/ha or in tCO2e/ha.

    Returns:
    float: Belowground biomass (Mg/ha) or Belowground carbon stock in tC/ha or in tCO2e/ha
    """
    return abg_b_or_c * ratio_below_ground_biomass_to_above_ground_biomass

def get_carbon_stock(data, scenario, carbon_time_id):
    return data.get("scenarios").get(scenario)[0][carbon_time_id]
  
#Define equation to calculate Extracted Log Emittions (tC/m3)
def calculate_ELE(D):
    """
    Calculate ELE
    
    Parameters:
    D = wood density (t/m3)
    
    Returns:
    float: ELE in tC/m3.
    """
    return (0.4924 * D) - 0.0158

#equation to calculate timber tree illegal
# def calculate_timber_tree_illegal(forest_area, illegal_logging_rate, ELE):
#     """
#     Calculate carbon in extracted timber from illegal logging in year n (tC)
#     
#     Parameters:
#     forest_area (float): Forest Area at project start (ha).
#     illegal_logging_rate (float): volume timber over bark extracted illegally each year (m3 ha-1 yr-1).
#     ELE (float): Extracted Log Emissions (t C m-3 extracted)
#     
#     Returns:
#     float: Carbon stock loss in tCO2e/yr.
#     """
#     forest_area = scenario_data.get("forest_area")
#     illegal_logging_rate = scenario_data.get("illegal_logging_rate")
#     
#     return forest_area * illegal_logging_rate * ELE

#equation to calculate logging damage factor (LDF)
def calculate_LDF(average_agbd_cstock):
    """
    Calculate LDF
    
    Parameters:
    average_agbd_cstock = total carbon stored in trees (above and belowground) in tC/ha.
    
    Returns:
    float: LDF in tC/ha.
    """
    return (-0.0039 * average_agbd_cstock) + 1.7817

#equation to calculate incidental damage
def calculate_incidental_damage(LDF, illegal_logging_rate):
    """
    Calculate incidental_damage
    
    Parameters:
    LDF in tC/ha.
    illegal_logging_rate in volume timber over bark extracted illegally each year (m3 ha-1 yr-1).
    
    Returns:
    float: incidental_damage in tC/m3.
    """
    return LDF * illegal_logging_rate
  
#equation to calculate historical deforestation
def calculate_historical_deforestation_rate(deforestation_rate_pre, deforestation_rate_post):
    """
    Calculate the historical deforestation rate.
    
    Parameters:
    deforestation_rate_pre (float): Deforestation rate before intervention (%/yr).
    deforestation_rate_post (float): Deforestation rate after intervention (%/yr).
    
    Returns:
    float: Historical deforestation rate (%/yr).
    """
    return deforestation_rate_pre - deforestation_rate_post

#equation to calculate deforested area
def calculate_area_deforested(historical_deforestation_rate, forest_area):
    """
    Calculate the area deforested.
    
    Parameters:
    historical_deforestation_rate (float): Historical deforestation rate (%/yr).
    forest_area (float): Forest area in ha.
    
    Returns:
    float: Area deforested in ha/yr.
    """
    return (historical_deforestation_rate / 100) * forest_area

#equation to calculate total carbon stock loss due to deforestation
def calculate_carbon_stock_loss(area_deforested, average_total_tco2e):
    """
    Calculate the carbon stock loss.
    
    Parameters:
    area_deforested (float): Area deforested in ha/yr.
    average_total_tco2e (float): Total (aboveground and belowground) biomass carbon stock in tCO2e/ha.
    
    Returns:
    float: Carbon stock loss in tCO2e/yr.
    """
    return area_deforested * average_total_tco2e

#equation to calculate annual baseline emissions
def calculate_annual_baseline_emissions(historical_deforestation_rate, carbon_stock_loss):
    """
    Calculate the annual baseline emissions.
    
    Parameters:
    historical_deforestation_rate (float): Historical deforestation rate (%/yr).
    carbon_stock_loss (float): Carbon stock loss in tCO2e/yr.
    
    Returns:
    float: Annual baseline emissions in tCO2e/yr.
    """
    return historical_deforestation_rate * carbon_stock_loss

def sanitize_filename(filename, max_length=45):
    # Replace spaces with underscores and remove any non-alphanumeric characters
    sanitized = re.sub(r'[^\w\s-]', '', filename.lower())
    sanitized = re.sub(r'[-\s]+', '_', sanitized)
    
    # Truncate to max_length characters
    return sanitized[:max_length]
