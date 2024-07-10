
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

# def get_scenario_parameter(data, scenario, param):
#     return data.get("scenarios").get(scenario)[0][param]
