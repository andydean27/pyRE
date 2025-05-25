from typing import Optional
import math
from pyRE.fluid.models import Oil, Water, Gas
from pyRE.fluid.collections import FluidCollection
from pyRE.rock.models import *
from pyRE.reservoir.correlations import *
from pyRE.reservoir.state import ReservoirState, Saturations

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
            fluids: Optional[FluidCollection],
            relative_permeability: Optional[dict[str, object]],
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
        self.rock = rock # Delegate rock-specific behavior to the Rock object
        self.fluids = fluids
        self.relative_permeability = relative_permeability or {}


        # Initialise initial conditions
        self.set_initial_conditions(
            pressure = kwargs.get("pressure", None),
            temperature = kwargs.get("temperature", None),
            saturations = kwargs.get("saturations", {})
        )
    
    def __repr__(self):
        # Return json representation of the reservoir
        return f"{self.__class__.__name__}({{depth: {self.depth}, thickness: {self.thickness}, net_to_gross: {self.net_to_gross}, area: {self.area}, rock: {self.rock}, fluids: {self.fluids}}})"    

    def set_initial_conditions(self,
                               pressure: float,
                               temperature: float,
                               saturations: dict[str, float]):
        """
        Set the initial conditions for the reservoir.
        
        Arguments:
        - pressure (float): Initial pressure [psia].
        - temperature (float): Initial temperature [Rankine].
        - saturations (dict): Dictionary of initial saturations for each fluid.
            - oil (float): Initial oil saturation.
            - gas (float): Initial gas saturation.
            - water (float): Initial water saturation.
        """
        self.initial_state = ReservoirState(
            pressure = pressure,
            temperature = temperature,
            saturations = Saturations(
                oil = saturations.get('oil', None),
                gas = saturations.get('gas', None),
                water = saturations.get('water', None)
            )
        )
    
    def verify_reservoir_setup(self):
        """
        Verify that the reservoir is set up correctly.
        This includes checking if the rock and fluids are set up correctly.
        """
        # Check if the rock is set
        if self.rock is None:
            raise ValueError("Rock properties must be set.")

        # Check if the fluids are set
        if not self.fluids:
            raise ValueError("Fluids must be set.")

        # Check if the initial conditions are set
        if self.initial_state is None:
            raise ValueError("Initial conditions must be set.")
        
        # Check saturation is set for each fluid present
        if self.fluids.oil is not None and self.initial_state.saturations.oil is None:
            raise ValueError("Initial oil saturation must be set if oil phase is present.")
        if self.fluids.gas is not None and self.initial_state.saturations.gas is None:
            raise ValueError("Initial gas saturation must be set if gas phase is present.")
        if self.fluids.water is not None and self.initial_state.saturations.water is None:
            raise ValueError("Initial water saturation must be set if water phase is present.")
        
        # Check if the initial saturations sum to 1
        if self.initial_state.saturations.sum() != 1.0:
            raise ValueError("Initial saturations must sum to 1.0.")
    

    def calculate_in_place_volumes(self, state: ReservoirState):
        """
        Calculate the in-place volumes of oil, gas, and water in the reservoir.
        
        Arguments:
        - state (ReservoirState): The state of the reservoir.
        
        Returns:
        - dict: A dictionary containing the in-place volumes of oil, gas, and water.
        """

        # Initialise in-place volumes
        oil_in_place = None
        gas_in_place = None
        water_in_place = None

        # Calculate in-place volumes based on the fluid present
        if self.fluids.oil is not None:
            oil_in_place = 0 # TODO

        if self.fluids.gas is not None:
            gas_in_place = self.pore_volume/1000000 * state.saturations.gas / self.fluids.gas.formation_volume_factor(state.pressure, state.temperature)
        
        if self.fluids.water is not None:
            water_in_place = self.pore_volume/5.615 * state.saturations.water / self.fluids.water.formation_volume_factor

        return oil_in_place, gas_in_place, water_in_place
        
    def calculate_initial_in_place_volumes(self):

        # Verify initial conditions are set
        if self.initial_state is None:
            raise ValueError("Initial conditions must be set.")

        # Calculate the in-place volumes of oil, gas, and water in the reservoir
        oil_in_place, gas_in_place, water_inplace = self.calculate_in_place_volumes(self.initial_state)

        # Store the initial in-place volumes
        self.initial_oil_in_place = oil_in_place
        self.initial_gas_in_place = gas_in_place
        self.initial_water_in_place = water_inplace
        


class CoalSeamGasReservoir(Reservoir):
    """
    Reservoir class for coal seam gas reservoirs.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Check if the rock is of type Coal
        if not isinstance(self.rock, Coal):
            raise ValueError("CoalSeamGasReservoir requires a CoalRock object.")
    
    def verify_reservoir_setup(self):
        super().verify_reservoir_setup()

        # Coal specific checks
        
        # Check if the rock is of type Coal
        if not isinstance(self.rock, Coal):
            raise ValueError("CoalSeamGasReservoir requires a CoalRock object.")
        
        # Check no oil properties are set
        if self.fluids.oil is not None:
            raise ValueError("CoalSeamGasReservoir does not support oil. Oil phase will be ignored.")
        if self.initial_state.saturations.oil is not None:
            raise ValueError("CoalSeamGasReservoir does not support oil. Remove oil saturation from initial conditions")
    
    def calculate_in_place_volumes(self, state):

        # Calculate adsorbed gas in place
        # 0.028316 converts g/cc to tonnes/cf
        # 43560 converts acres to cf
        # 1000000 converts cf to mmcf
        adsorbed_gas_in_place = 0.028316 * 43560 / 1000000 * self.area * self.thickness * self.net_to_gross * self.rock.density * self.rock.gas_content * (1 - self.rock.ash_fraction - self.rock.moisture_fraction)

        oil_in_place, gas_in_place, water_in_place = super().calculate_in_place_volumes(state)
        return oil_in_place, gas_in_place, water_in_place, adsorbed_gas_in_place
    
    def calculate_initial_in_place_volumes(self):

        # Verify initial conditions are set
        if self.initial_state is None:
            raise ValueError("Initial conditions must be set.")

        # Calculate the in-place volumes of oil, gas, and water in the reservoir
        oil_in_place, gas_in_place, water_inplace, adsorbed_gas_in_place = self.calculate_in_place_volumes(self.initial_state)

        # Store the initial in-place volumes
        self.initial_oil_in_place = oil_in_place
        self.initial_gas_in_place = gas_in_place
        self.initial_water_in_place = water_inplace
        self.initial_adsorbed_gas_in_place = adsorbed_gas_in_place
    
    def sorption_compressibility(self, pressure: float = None, temperature: float = None):
        """
        Calculate the sorption compressibility for coal seam gas reservoir.
        """
        # Use initial conditions if no inputs are provided
        pressure = pressure or self.initial_state.pressure
        temperature = temperature or self.initial_state.temperature

        # Calculate sorption compressibility using the Langmuir isotherm
                # 0.0283168 is a conversion from g/cc to tonne/scf
        sorption_compressibility = 0.028316*(
            (self.fluids.gas.formation_volume_factor(pressure, temperature)*self.rock.langmuir_volume*(1-self.rock.ash_fraction-self.rock.moisture_fraction)*self.rock.density*self.rock.langmuir_pressure)/
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
        pressure = pressure or self.initial_state.pressure
        temperature = temperature or self.initial_state.temperature
        water_saturation = water_saturation or self.initial_state.saturations.water

        # Calculate Z* using Modified King's method
        z_star = (
            self.fluids.gas.z(pressure, temperature) / 
                ((self.rock.density * self.rock.langmuir_volume * (1 - self.rock.ash_fraction - self.rock.moisture_fraction) * 14.7 * temperature * self.fluids.gas.z(pressure, temperature)) / 
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
                Options are "modified_king", "jensen_smith", clarkson_mcgovern.

        """
        from sklearn.linear_model import LinearRegression
        valid_methods = ["modified_king", "jensen_smith", "clarkson_mcgovern"]

        # Check if the required columns are present in the DataFrame
        required_columns = ["pressure", "cumulative_gas_production", "cumulative_water_production"]
        for col in required_columns:
            if col not in data.columns:
                raise ValueError(f"Missing required column: {col}")

        # Calculate the material balance using the specified method
        if method == "modified_king":
            result = self._modified_king_material_balance(data)
        elif method == "jensen_smith":
            result = self._jensen_smith_material_balance(data)
        elif method == "clarkson_mcgovern":
            result = self._clarkson_mcgovern_material_balance(data)
        else:
            raise ValueError(f"Invalid method: {method}. Valid options are: {', '.join(valid_methods)}")
        
        # Perform linear regression on the material balance data
        y_identifier = {
            "modified_king": "p/z*",
            "jensen_smith": "p/(p+p_L)",
            "clarkson_mcgovern": "C&M MBE"
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

    
    def _modified_king_water_saturation(self, cumulative_water_production: float):
        """
        Calculate the water saturation using Modified King's method. (largely ignores water and formation compressibilityies)
        
        Arguments:
        - cumulative_water_production (float): Cumulative water production [bbl].
        
        Returns:
        - water_saturation (float): Water saturation.
        """
        # Calculate water required to produce to reach critical desorption pressure
        cumulative_water_to_desorption = self.pore_volume/5.615 * self.rock.compressibility * max(0, (self.initial_state.pressure - self.rock.critical_desorption_pressure))

        # Calculate water saturation
        water_saturation = (
            self.initial_state.saturations.water - 
            self.fluids.water.formation_volume_factor/(self.pore_volume/5.615)*max(0, (cumulative_water_production - cumulative_water_to_desorption))
        )
        
        return water_saturation
    
    def _modified_king_material_balance(self, data: pd.DataFrame):
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
        result["water_saturation"] = result.apply(lambda row: self._modified_king_water_saturation(row["cumulative_water_production"]), axis=1)
        # Calculate Z*
        result["z*"] = result.apply(lambda row: self.z_star(pressure = row["pressure"],water_saturation = row["water_saturation"]), axis=1)
        # Calculate P/z*
        result["p/z*"] = result["pressure"] / result["z*"]
        
        return result
    
    def _jensen_smith_material_balance(self, data: pd.DataFrame):
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

    def _clarkson_mcgovern_material_balance(self, data: pd.DataFrame):
        """
        Calculate gas material balance using Clarkson and McGovern method.

        Arguments:
        - data (pd.DataFrame): DataFrame containing the observed or simulated data.
            - pressure (float): Pressure [psia].
            - cumulative_gas_production (float): Cumulative gas production [MMscf].
            - cumulative_water_production (float): Cumulative water production [bbl].
        """

        result = data.copy()

        # Calculate C&M
        result["C&M MBE"] = result.apply(lambda row: (row["pressure"] / (row["pressure"] + self.rock.langmuir_pressure)) + 
                                         (32.037 * self.rock.porosity * (1-self.initial_state.saturations.water)/self.rock.langmuir_volume/self.fluids.gas.formation_volume_factor(row["pressure"], self.initial_state.temperature)/self.rock.density), axis=1)

        return result