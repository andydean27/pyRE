import pytest
from pyRE.reservoir.models import CoalSeamGasReservoir
from pyRE.rock.models import Coal
from pyRE.fluid.models import Oil, Gas, Water

def test_coal_seam_gas_reservoir_initialization():
    # Create a CoalRock object
    rock = Coal(porosity=0.1, permeability=5, compressibility=1e-5, gas_content=500, density = 1.75, langmuir_pressure = 1000, langmuir_volume = 200)

    # Create fluid objects
    oil = Oil()
    gas = Gas()
    water = Water()

    # Initialize the CoalSeamGasReservoir
    reservoir = CoalSeamGasReservoir(
        depth=1000,
        thickness=50,
        net_to_gross=0.8,
        area=200,
        rock=rock,
        oil=oil,
        gas=gas,
        water=water,
        skin = 0
    )

    # Assertions to verify initialization
    assert reservoir.depth == 1000
    assert reservoir.thickness == 50
    assert reservoir.net_to_gross == 0.8
    assert reservoir.area == 200
    assert reservoir.rock == rock
    assert reservoir.oil == oil
    assert reservoir.gas == gas
    assert reservoir.water == water
    assert reservoir.skin == 0