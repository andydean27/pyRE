from typing import Optional
from pyRE.fluid.models import Oil, Water, Gas
from pyRE.rock.models import *

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
                pressure: float,
                temperature: float,
                rock: Rock,
                oil_saturation: float = 1.0,
                gas_saturation: float = 0.0,
                water_saturation: float = 0.0,
                oil: Optional[Oil] = None,
                gas: Optional[Gas] = None,
                water: Optional[Water] = None,
                skin: float = 0
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
        self.pressure = pressure
        self.temperature = temperature
        self.rock = rock  # Delegate rock-specific behavior to the Rock object
        self.oil = oil or Oil()  # Default to a new Oil object if not provided
        self.gas = gas or Gas()  # Default to a new Gas object if not provided
        self.water = water or Water()  # Default to a new Water object if not provided

        # Validate saturations
        if oil_saturation + gas_saturation + water_saturation != 1.0:
            raise ValueError("Total saturation must equal 1.0")
        self.oil_saturation = oil_saturation
        self.gas_saturation = gas_saturation
        self.water_saturation = water_saturation
        
    def calculate_in_place_volumes(self):
        oil_in_place = 0 # TODO

        gas_in_place = 0 # TODO

        water_inplace = self.area * self.thickness * self.net_to_gross * self.rock.porosity * self.water_saturation / self.water.formation_volume_factor / 5.625 # Convert to bbl

        return {
                "oil_in_place": oil_in_place,
                "gas_in_place": gas_in_place,
                "water_in_place": water_inplace
                }
    
    def calculate_total_compressibility(self):
        pass

class ConventionalReservoir(Reservoir):
    """
    Reservoir class for conventional reservoirs.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class CoalSeamGasReservoir(Reservoir):
    """
    Reservoir class for coal seam gas reservoirs.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Check if the rock is of type CoalRock
        if not isinstance(self.rock, CoalRock):
            raise ValueError("CoalSeamGasReservoir requires a CoalRock object.")
        
        # Validate saturations
        if self.oil_saturation != 0.0:
            raise ValueError("Oil saturation must be 0.0 for CoalSeamGasReservoir.")

    def calculate_adsorbed_gas_in_place(self):
        pass