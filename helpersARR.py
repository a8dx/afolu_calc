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
    scenario_data = data["scenarios"][scenario][0]
    aoi_subregions = scenario_data["aoi_subregions"]
    
    carbon_stocks = []
    
    for subregion in aoi_subregions:
        aoi_id = subregion["aoi_id"]
        area = subregion["area"]
        
        # Get the carbon stock value for this subregion
        # If the carbon_time_id is not specific to subregions, we'll use the scenario-level value
        carbon_stock = subregion.get(carbon_time_id, scenario_data.get(carbon_time_id))
        
        if carbon_stock is None:
            print(f"Warning: Carbon stock data '{carbon_time_id}' not found for subregion {aoi_id}")
            carbon_stock = 0
        
        carbon_stocks.append({
            "aoi_id": aoi_id,
            "area": area,
            "carbon_stock": float(carbon_stock)
        })
        
        print(f"Subregion {aoi_id}: Carbon stock ({carbon_time_id}) = {carbon_stock}")
    
    return carbon_stocks

def sanitize_filename(filename, max_length=45):
    # Replace spaces with underscores and remove any non-alphanumeric characters
    sanitized = re.sub(r'[^\w\s-]', '', filename.lower())
    sanitized = re.sub(r'[-\s]+', '_', sanitized)
    
    # Truncate to max_length characters
    return sanitized[:max_length]
