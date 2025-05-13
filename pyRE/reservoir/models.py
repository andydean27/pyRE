from typing import Optional
import math
from pyRE.fluid.models import Oil, Water, Gas
from pyRE.rock.models import *
from pyRE.reservoir.correlations import *

import pandas as pd
import numpy as np

class Reservoir:
    """
    Base class for reservoirs.
    """
    def __init__(
            self,
            depth: float,
            thickness: float,
            net_to_gross: float,
            area: float,
            rock: Optional[Rock],
            oil: Optional[Oil] = None,
            gas: Optional[Gas] = None,
            water: Optional[Water] = None,
            water_oil_relative_permeability = None,
            gas_oil_relative_permeability = None,
            water_gas_relative_permeability = None,
            skin: float = 0,
            **kwargs
    ):
        """
        Initialize a reservoir.

        Arguments:
        - depth (float): Depth [ft].
        - thickness (float): Thickness [ft].
        - net_to_gross (float): Net-to-gross ratio.
        - area (float): Area [acres].
        - pressure (float): Pressure [psia].
        - temperature (float): Temperature [Rankine].
        - rock (Rock): Rock properties object.
        - oil (Oil, optional): Oil object.
        - gas (Gas, optional): Gas object.
        - water (Water, optional): Water object.
        """
        self.depth = depth
        self.thickness = thickness
        self.net_to_gross = net_to_gross
        self.area = area
        self.bulk_volume = 43560 * area * thickness * net_to_gross  # Bulk volume of the reservoir [ft3]
        self.pore_volume = self.bulk_volume * rock.porosity  # Pore volume of the reservoir [ft3]
        self.rock = rock or Conventional() # Delegate rock-specific behavior to the Rock object
        self.oil = oil or Oil()  # Default to a new Oil object if not provided
        self.gas = gas or Gas()  # Default to a new Gas object if not provided
        self.water = water or Water()  # Default to a new Water object if not provided
        self.water_oil_relative_permeability = water_oil_relative_permeability
        self.gas_oil_relative_permeability = gas_oil_relative_permeability
        self.water_gas_relative_permeability = water_gas_relative_permeability
        self.skin = skin


        # Initialise initial conditions
        self.initial_pressure = None
        self.initial_temperature = None
        self.initial_oil_saturation = 0
        self.initial_gas_saturation = 0
        self.initial_water_saturation = 0

    def set_initial_conditions(self, **kwargs):
        initial_conditions = [
            "initial_pressure",
            "initial_temperature",
            "initial_oil_saturation",
            "initial_gas_saturation",
            "initial_water_saturation"
        ]

        for key, value in kwargs.items():
            if key in initial_conditions:
                setattr(self, key, value)
            else:
                raise ValueError(f"Invalid initial condition: {key}. Valid options are: {', '.join(initial_conditions)}")
        
        # Validate saturations
        if self.initial_oil_saturation + self.initial_gas_saturation + self.initial_water_saturation != 1.0:
            raise ValueError("Saturations must sum to 1.0.")
        

        
    def calculate_initial_in_place_volumes(self):

        # Verify initial conditions are set

        oil_in_place = 0 # TODO

        gas_in_place = 0.043560 * self.area * self.thickness * self.net_to_gross * self.rock.porosity * self.initial_gas_saturation / self.gas.formation_volume_factor(self.initial_pressure, self.initial_temperature)

        water_inplace = self.pore_volume/5.615 * self.initial_water_saturation / self.water.formation_volume_factor

        self.initial_oil_in_place = oil_in_place
        self.initial_gas_in_place = gas_in_place
        self.initial_water_in_place = water_inplace
        
        return {
                "oil_in_place": oil_in_place,
                "gas_in_place": gas_in_place,
                "water_in_place": water_inplace
                }


class CoalSeamGasReservoir(Reservoir):
    """
    Reservoir class for coal seam gas reservoirs.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Check if the rock is of type Coal
        if not isinstance(self.rock, Coal):
            raise ValueError("CoalSeamGasReservoir requires a CoalRock object.")
        
        # Validate saturations
        if self.initial_oil_saturation != 0.0:
            raise ValueError("Oil saturation must be 0.0 for CoalSeamGasReservoir.")

    def calculate_initial_in_place_volumes(self):

        adsorbed_gas_in_place = 0.0013597 * self.area * self.thickness * self.net_to_gross * self.rock.density * self.rock.gas_content * (1 - self.rock.ash_fraction - self.rock.moisture_fraction)
        self.initial_adsorbed_gas_in_place = adsorbed_gas_in_place

        return {
            **super().calculate_initial_in_place_volumes(),
            "adsorbed_gas_in_place": adsorbed_gas_in_place
        }
    
    def sorption_compressibility(self, pressure: float = None, temperature: float = None):
        """
        Calculate the sorption compressibility for coal seam gas reservoir.
        """
        # Use initial conditions if no inputs are provided
        pressure = pressure or self.initial_pressure
        temperature = temperature or self.initial_temperature

        # Calculate sorption compressibility using the Langmuir isotherm
                # 0.02787 is a conversion from g/cc to ton/scf
        sorption_compressibility = 0.02787*(
            (self.gas.formation_volume_factor(pressure, temperature)*self.rock.langmuir_volume*(1-self.rock.ash_fraction-self.rock.moisture_fraction)*self.rock.density*self.rock.langmuir_pressure)/
            (self.rock.porosity*(pressure + self.rock.langmuir_pressure)**2)
        )
        
        return sorption_compressibility

    def z_star(self, pressure: float = None, temperature: float = None, water_saturation: float = None):
        """
        Calculate Modified King's Z* for coal seam gas reservoirs.
        If no inputs are provided use initial conditions.
        
        Arguments:
        - pressure (float): Pressure [psia].
        - temperature (float): Temperature [Rankine].
        """
        # Use initial conditions if no inputs are provided
        pressure = pressure or self.initial_pressure
        temperature = temperature or self.initial_temperature
        water_saturation = water_saturation or self.initial_water_saturation

        # Calculate Z* using Modified King's method
        z_star = (
            self.gas.z(pressure, temperature) / 
                ((self.rock.density * self.rock.langmuir_volume * (1 - self.rock.ash_fraction - self.rock.moisture_fraction) * 14.7 * temperature * self.gas.z(pressure, temperature)) / 
                 ((32.037*self.rock.porosity*536.7*(pressure+self.rock.langmuir_pressure)) + (1 - water_saturation)))
        )
        
        return z_star
    
    ##############################
    # Material balance methods
    ##############################


    def gas_material_balance(self, data: pd.DataFrame, method: str = "modified_king"):
        """
        Calculate the material balance data and plot (optional) for a set of given observed or simulated data
        
        Arguments:
        - data (pd.DataFrame): DataFrame containing the observed or simulated data.
            - pressure (float): Pressure [psia].
            - cumulative_gas_production (float): Cumulative gas production [MMscf].
            - cumulative_water_production (float): Cumulative water production [bbl].
        
        - method (str): Method to use for the material balance calculation. 
                Options are "modified_king", "jensen_smith".

        """
        from sklearn.linear_model import LinearRegression
        valid_methods = ["modified_king", "jensen_smith"]

        # Check if the required columns are present in the DataFrame
        required_columns = ["pressure", "cumulative_gas_production", "cumulative_water_production"]
        for col in required_columns:
            if col not in data.columns:
                raise ValueError(f"Missing required column: {col}")

        # Calculate the material balance using the specified method
        if method == "modified_king":
            result = self.modified_king_material_balance(data)
        elif method == "jensen_smith":
            result = self.jensen_smith_material_balance(data)
        else:
            raise ValueError(f"Invalid method: {method}. Valid options are: {', '.join(valid_methods)}")
        
        # Perform linear regression on the material balance data
        y_identifier = {
            "modified_king": "p/z*",
            "jensen_smith": "p/(p+p_L)"
        }
        model = LinearRegression()
        model.fit(result['cumulative_gas_production'].values.reshape(-1, 1), result[y_identifier[method]].values.reshape(-1, 1))
        # Get x intersect
        x_intercept = -model.intercept_[0] / model.coef_[0][0]
        # Calculate gas in place difference
        gas_in_place_variance = self.initial_gas_in_place - x_intercept
        
        return result, {
            "model": model,
            "gas_in_place": x_intercept
        }

    
    def modified_king_water_saturation(self, cumulative_water_production: float):
        """
        Calculate the water saturation using Modified King's method. (largely ignores water and formation compressibilityies)
        
        Arguments:
        - cumulative_water_production (float): Cumulative water production [bbl].
        
        Returns:
        - water_saturation (float): Water saturation.
        """
        # Calculate water required to produce to reach critical desorption pressure
        cumulative_water_to_desorption = self.pore_volume/5.615 * self.rock.compressibility * max(0, (self.initial_pressure - self.rock.critical_desorption_pressure))

        # Calculate water saturation
        water_saturation = (
            self.initial_water_saturation - 
            self.water.formation_volume_factor/(self.pore_volume/5.615)*max(0, (cumulative_water_production - cumulative_water_to_desorption))
        )
        
        return water_saturation
    
    def modified_king_material_balance(self, data: pd.DataFrame):
        """
        Calculate gas material balance using Modified King's method.

        Arguments:
        - data (pd.DataFrame): DataFrame containing the observed or simulated data.
            - pressure (float): Pressure [psia].
            - cumulative_gas_production (float): Cumulative gas production [MMscf].
            - cumulative_water_production (float): Cumulative water production [bbl].
        """
        
        result = data.copy()
        # Calculate water saturation
        result["water_saturation"] = result.apply(lambda row: self.modified_king_water_saturation(row["cumulative_water_production"]), axis=1)
        # Calculate Z*
        result["z*"] = result.apply(lambda row: self.z_star(pressure = row["pressure"],water_saturation = row["water_saturation"]), axis=1)
        # Calculate P/z*
        result["p/z*"] = result["pressure"] / result["z*"]
        
        return result
    
    def jensen_smith_material_balance(self, data: pd.DataFrame):
        """
        Calculate gas material balance using Jensen and Smith method.

        Arguments:
        - data (pd.DataFrame): DataFrame containing the observed or simulated data.
            - pressure (float): Pressure [psia].
            - cumulative_gas_production (float): Cumulative gas production [MMscf].
            - cumulative_water_production (float): Cumulative water production [bbl].
        """
        result = data.copy()
        # Calculate p/(p+pl)
        result["p/(p+p_L)"] = result.apply(lambda row: row["pressure"] / (row["pressure"] + self.rock.langmuir_pressure), axis=1)

        
        return result