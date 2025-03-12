import pytest
from pyRE.fluids.models import *

#------------------------------------------------------------
#   Water Tests
#------------------------------------------------------------

def test_water_initialization():
    w = water(Bw=1.0, mu=1.0, rho=1.0)
    assert w.Bw == 1.0
    assert w.mu == 1.0
    assert w.rho == 1.0

#------------------------------------------------------------
#   Gas Tests
#------------------------------------------------------------

def mock_z_correlation(P, T, gas_obj):
    # Mock z-factor correlation function
    return 1.0

def mock_mu_correlation(P, T, gas_obj):
    # Mock viscosity correlation function
    return 0.01

def test_gas_initialization():
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=mock_z_correlation, mu_correlation=mock_mu_correlation)
    assert g.Pc == 1000
    assert g.Tc == 500
    assert g.M == 20
    assert g.z_correlation == mock_z_correlation
    assert g.mu_correlation == mock_mu_correlation

def test_gas_z():
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=mock_z_correlation, mu_correlation=mock_mu_correlation)
    assert g.z(1000, 500) == 1.0

def test_gas_mu():
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=mock_z_correlation, mu_correlation=mock_mu_correlation)
    assert g.mu(1000, 500) == 0.01

def test_gas_Bg():
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=mock_z_correlation, mu_correlation=mock_mu_correlation)
    assert g.Bg(1000, 500) == pytest.approx(0.0283 * 500 / 1000, rel=1e-3)

def test_gas_rho():
    g = gas(Pc=1000, Tc=500, M=20, z_correlation=mock_z_correlation, mu_correlation=mock_mu_correlation)
    assert g.rho(1000, 500) == pytest.approx(144 * 1000 * 20 / 1.0 / 1545.349 / 500, rel=1e-3)