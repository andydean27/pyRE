

import pytest
from pyRE.fluid.models import gas
from pyRE.fluid.correlations import STANDINGKATZ, LEEGONZALEZEAKIN

def test_STANDINGKATZ():
    # Create a mock gas object
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=STANDINGKATZ, mu_correlation=LEEGONZALEZEAKIN)
    
    # Test the STANDINGKATZ function
    P = 100
    T = 500
    expected_z = STANDINGKATZ(P, T, g)
    assert expected_z == pytest.approx(0.9676, rel=1e-2) 

def test_LEEGONZALEZEAKIN():
    # Create a mock gas object
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=STANDINGKATZ, mu_correlation=LEEGONZALEZEAKIN)
    
    # Test the LEEGONZALEZEAKIN function
    P = 100
    T = 500
    expected_mu = LEEGONZALEZEAKIN(P, T, g)
    assert expected_mu == pytest.approx(0.0101, rel=1e-2) 