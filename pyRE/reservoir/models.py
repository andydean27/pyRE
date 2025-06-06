from typing import Optional
import math
from pyRE.fluid.models import Oil, Water, Gas
from pyRE.fluid.collections import FluidCollection
from pyRE.rock.models import *
from pyRE.reservoir.state import ReservoirState, Saturations
from pyRE.relative_permeability.models import RelativePermeability

import pandas as pd
import numpy as np


shape_factors = {
    'circular': 31.62,
    'hexagonal': 31.6,
    'square': 30.8828,
    'equilateral_triangle': 27.6,
}


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
            rock: Rock,
            fluids: FluidCollection,
            relative_permeability: Optional[RelativePermeability],
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
        - fluids (FluidCollection): Collection of fluid models (Oil, Gas, Water).
        - relative_permeability (RelativePermeability, optional): Relative permeability model.

        """
        self.depth = depth
        self.thickness = thickness
        self.net_to_gross = net_to_gross
        self.area = area
        self.bulk_volume = 43560 * area * thickness * net_to_gross  # Bulk volume of the reservoir [ft3]
        self.pore_volume = self.bulk_volume * rock.porosity  # Pore volume of the reservoir [ft3]
        self.rock = rock # Delegate rock-specific behavior to the Rock object
        self.fluids = fluids
        self.relative_permeability = relative_permeability


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
        
    def _validate_realtive_permeability(self):
        """
        Validate that the relative permeability model is set up correctly.
        This includes checking if the relative permeability model is set, if it is of the correct type and if the curves set match the fluids present.
        """
        if self.relative_permeability is None:
            raise ValueError("Relative permeability model must be set for any flow simulations.")
        
        if not isinstance(self.relative_permeability, RelativePermeability):
            raise TypeError(f"Expected RelativePermeability model, got {type(self.relative_permeability)}")
        
        # Check if the relative permeability curves are set for each fluid present
        if self.fluids.oil is not None and self.relative_permeability.oil_curve is None:
            raise ValueError("Relative permeability curve must be set for oil phase if oil is present")
        if self.fluids.oil is None and self.relative_permeability.oil_curve is not None:
            raise ValueError("Relative permeability curve for oil phase is set, but oil phase is not present in the reservoir")
        if self.fluids.gas is not None and self.relative_permeability.gas_curve is None:
            raise ValueError("Relative permeability curve must be set for gas phase if gas is present")
        if self.fluids.gas is None and self.relative_permeability.gas_curve is not None:
            raise ValueError("Relative permeability curve for gas phase is set, but gas phase is not present in the reservoir")
        if self.fluids.water is not None and self.relative_permeability.water_curve is None:
            raise ValueError("Relative permeability curve must be set for water phase if water is present")
        if self.fluids.water is None and self.relative_permeability.water_curve is not None:
            raise ValueError("Relative permeability curve for water phase is set, but water phase is not present in the reservoir")
    
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
            water_in_place = self.pore_volume/5.615 * state.saturations.water / self.fluids.water.formation_volume_factor(state.pressure, state.temperature)

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
        

    # Vertical well inflow equations
    def _pseudo_steady_flow_rate(
            self,
            state: ReservoirState,
            fluid: str,
            well_flowing_pressure: float,
            wellbore_radius: float,
            skin: float,
            shape_factor: str = 'circular',
            pressure_squared_approximation: bool = True) -> float:
        """
        Calculate the pseudo-steady flow rate for a vertical well.
        
        Arguments:
        - state: ReservoirState object containing current state of the reservoir (pressure, temperature, saturation).
        - fluid: fluid phase for which the rate is calculate ('oil', 'gas', or 'water').
        - well_flowing_pressure: Flowing pressure at the wellbore in psi.
        - wellbore_radius: Radius of the wellbore in feet.
        - skin: Skin factor of the well.
        - shape_factor: Shape factor for the wellbore, default is 31.62 for circular wells.
            - Can also be a string ('circular', 'hexagonal', 'square', 'equilateral_triangle') or a float value.
            - https://onepetro.org/spe/general-information/1577/Fluid-flow-through-permeable-media?searchresult=1
        - pressure_squared_approximation: If True, uses the pressure squared approximation for gas pseudo pressure.
        
        Returns:
        - Flow rate in STB/day for oil, SCF/day for gas, or BWPD for water.
        """


        shape_factor = shape_factors.get(shape_factor.lower(), 31.62)
        if shape_factor is None:
            raise ValueError(f"Invalid shape factor: {shape_factor}")
 
        
        # Check if flowing pressure is below reservoir pressure
        if well_flowing_pressure > state.pressure:
            return 0.0


        # Get the fluid model
        if fluid == 'oil':
            if self.fluids.oil is None:
                raise ValueError("Oil model is not defined in the reservoir.")
            fluid = self.fluids.oil
        elif fluid == 'gas':
            if self.fluids.gas is None:
                raise ValueError("Gas model is not defined in the reservoir.")
            fluid = self.fluids.gas
        elif fluid == 'water':
            if self.fluids.water is None:
                raise ValueError("Water model is not defined in the reservoir.")
            fluid = self.fluids.water
        else:
            raise ValueError(f"Invalid fluid type: {fluid}. Must be 'oil', 'gas', or 'water'.")

        # Calculate relative permeabilities
        self.relative_permeability.calculate(state.saturations)
        if fluid == self.fluids.oil:
            relative_permeability = self.relative_permeability.oil
        elif fluid == self.fluids.gas:
            relative_permeability = self.relative_permeability.gas
        elif fluid == self.fluids.water:
            relative_permeability = self.relative_permeability.water
        else:
            raise ValueError(f"Fluid {fluid} not found in the reservoir's relative permeability model.")
        
        # Calculate for gas
        if isinstance(fluid, Gas):
            if pressure_squared_approximation:
                return (
                    ((state.pressure**2 - well_flowing_pressure**2) * self.rock.permeability * relative_permeability * self.thickness * self.net_to_gross) /
                    (1422 * state.temperature * fluid.z(state.pressure, state.temperature) * fluid.viscosity(state.pressure, state.temperature)) / 
                    (0.5 * math.log(10.06 * self.area * 43560 / (shape_factor * wellbore_radius**2)) - 0.75 + skin)
                )
            else:
                raise NotImplementedError("pseudo pressure not implemented yet")
            
        # Calculate for oil or water
        return (
            ((state.pressure - well_flowing_pressure) * self.rock.permeability * relative_permeability * self.thickness * self.net_to_gross) /
            (141.2 * fluid.formation_volume_factor(state.pressure, state.temperature) * fluid.viscosity(state.pressure, state.temperature)) / 
            (0.5 * math.log(10.06 * self.area * 43560 / (shape_factor * wellbore_radius**2)) - 0.75 + skin)
        )

    def inflow_performance_relationship(
            self,
            state: ReservoirState = None,
            flow_behaviour: str = "pseudo_steady",
            **kwargs
            ):
        """
        Calculate the inflow performance relationship for the reservoir.
        
        Arguments:
        - state (ReservoirState): The state of the reservoir. If not provided, uses the initial state.
        - flow_behaviour (str): The flow behaviour to use for the calculation.
            Supported: 'pseudo_steady'.
        - kwargs: Additional keyword arguments for the inflow equation."
            - wellbore_radius (float): Radius of the wellbore [ft]. Default is 0.5 ft.
            - skin (float): Skin factor of the well. Default is 0.
            - shape_factor (str | float): Shape factor for the wellbore. Default is 'circular'.
                Can also be a string ('circular', 'hexagonal', 'square', 'equilateral_triangle') or a float value.
            - pressure_squared_approximation (bool): If True, uses the pressure squared approximation for gas pseudo pressure. Default is True.
        
        Returns:
        - pd.DataFrame: DataFrame containing the inflow performance relationship.
            - pressure (float): Flowing pressure at the wellbore [psia].
            - oil_flow_rate (float): Flow rate of oil [STB/day].
            - gas_flow_rate (float): Flow rate of gas [SCF/day].
            - water_flow_rate (float): Flow rate of water [BWPD].
        """
        # Verify reservoir setup
        self.verify_reservoir_setup()

        # If state is not provided, use initial state
        state = state or self.initial_state

        # Get the inflow equation from flow behaviour
        if flow_behaviour == "pseudo_steady":
            inflow_equation = self._pseudo_steady_flow_rate
        else:
            raise ValueError(f"Invalid flow behaviour: {flow_behaviour}. Supported: 'pseudo_steady'.")
        
        # Generate flowing pressure range
        n = kwargs.get("n", 20)
        pressure_range = np.linspace(
            0, 
            state.pressure, 
            num=n
        )

        # Calculate flow rates for each pressure in the range for each phase present in the reservoir
        flow_rates = {}
        for fluid in ['oil', 'gas', 'water']:
            if getattr(self.fluids, fluid) is not None:
                flow_rates[fluid] = [
                    inflow_equation(
                        state=state,
                        fluid=fluid,
                        well_flowing_pressure=pressure,
                        wellbore_radius=kwargs.get("wellbore_radius", 0.5),
                        skin=kwargs.get("skin", 0),
                        shape_factor=kwargs.get("shape_factor", "circular"),
                        pressure_squared_approximation=kwargs.get("pressure_squared_approximation", True)
                    ) for pressure in pressure_range
                ]

        # Create a DataFrame to hold the results
        results = pd.DataFrame({
            "pressure": pressure_range,
            **{f"{fluid}_rate": flow_rates[fluid] for fluid in flow_rates}
        })
        # Set index to pressure
        results.set_index("pressure", inplace = True)

        return results

            

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
    

    # Material balance methods

    def gas_material_balance(self, data: pd.DataFrame, method: str = "modified_king"):
        """
        Calculate the material balance data and plot (optional) for a set of given observed or simulated data
        
        Arguments:
        - data (pd.DataFrame): DataFrame containing the observed or simulated data.
            - pressure (float): Pressure [psia].
            - temperature (float, optional): Temperature [Rankine]. If not provided, uses initial temperature.
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
            
        # Fill temperature with initial temperature if not provided
        if "temperature" not in data.columns:
            data["temperature"] = self.initial_state.temperature

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
    
    def _modified_king_water_saturation(self, cumulative_water_production: float, pressure: float = None, temperature: float = None):
        """
        Calculate the water saturation using Modified King's method. (largely ignores water and formation compressibilityies)
        
        Arguments:
        - cumulative_water_production (float): Cumulative water production [bbl].
        
        Returns:
        - water_saturation (float): Water saturation.
        """
        # Use initial conditions if no inputs are provided
        pressure = pressure or self.initial_state.pressure
        temperature = temperature or self.initial_state.temperature

        # Calculate water required to produce to reach critical desorption pressure
        cumulative_water_to_desorption = self.pore_volume/5.615 * self.rock.compressibility * max(0, (self.initial_state.pressure - self.rock.critical_desorption_pressure))

        # Calculate water saturation
        water_saturation = (
            self.initial_state.saturations.water - 
            self.fluids.water.formation_volume_factor(pressure, temperature)/(self.pore_volume/5.615)*max(0, (cumulative_water_production - cumulative_water_to_desorption))
        )
        
        return water_saturation
    
    def _modified_king_material_balance(self, data: pd.DataFrame):
        """
        Calculate gas material balance using Modified King's method.

        Arguments:
        - data (pd.DataFrame): DataFrame containing the observed or simulated data.
            - pressure (float): Pressure [psia].
            - temperature (float): Temperature [Rankine].
            - cumulative_gas_production (float): Cumulative gas production [MMscf].
            - cumulative_water_production (float): Cumulative water production [bbl].
        """
        
        result = data.copy()
        # Calculate water saturation
        result["water_saturation"] = result.apply(lambda row: self._modified_king_water_saturation(row["cumulative_water_production"], row["pressure"] or None, row["temperature"] or None), axis=1)
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
            - temperature (float): Temperature [Rankine].
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
            - temperature (float): Temperature [Rankine].
            - cumulative_gas_production (float): Cumulative gas production [MMscf].
            - cumulative_water_production (float): Cumulative water production [bbl].
        """

        result = data.copy()

        # Calculate C&M
        result["C&M MBE"] = result.apply(lambda row: (row["pressure"] / (row["pressure"] + self.rock.langmuir_pressure)) + 
                                         (32.037 * self.rock.porosity * (1-self.initial_state.saturations.water)/self.rock.langmuir_volume/self.fluids.gas.formation_volume_factor(row["pressure"], self.initial_state.temperature)/self.rock.density), axis=1)

        return result

