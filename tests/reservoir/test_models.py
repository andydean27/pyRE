import pytest
from pyRE.reservoir.models import CoalSeamGasReservoir
from pyRE.rock.models import CoalRock
from pyRE.fluid.models import Oil, Gas, Water

def test_coal_seam_gas_reservoir_initialization():
    # Create a CoalRock object
    rock = CoalRock(porosity=0.1, permeability=5, compressibility=1e-5, gas_content=500)

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
        pressure=3000,
        temperature=520,
        rock=rock,
        oil=oil,
        gas=gas,
        water=water
    )

    # Assertions to verify initialization
    assert reservoir.depth == 1000
    assert reservoir.thickness == 50
    assert reservoir.net_to_gross == 0.8
    assert reservoir.area == 200
    assert reservoir.pressure == 3000
    assert reservoir.temperature == 520
    assert reservoir.rock == rock
    assert reservoir.oil == oil
    assert reservoir.gas == gas
    assert reservoir.water == water