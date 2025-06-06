

import pytest
from pyRE.fluid.models import Gas
from pyRE.fluid.correlations import StandingKatz, LeeGonzalezEakin

def test_STANDINGKATZ():
    # Create a mock gas object
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=StandingKatz, viscosity_correlation=LeeGonzalezEakin)
    
    # Test the STANDINGKATZ function
    P = 100
    T = 500
    expected_z = StandingKatz(P, T, g)
    assert expected_z == pytest.approx(0.9676, rel=1e-2) 

def test_LEEGONZALEZEAKIN():
    # Create a mock gas object
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=StandingKatz, viscosity_correlation=LeeGonzalezEakin)
    
    # Test the LEEGONZALEZEAKIN function
    P = 100
    T = 500
    expected_mu = LeeGonzalezEakin(P, T, g)
    assert expected_mu == pytest.approx(0.0101, rel=1e-2) 