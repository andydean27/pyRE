import pytest
from pyRE.fluid.models import *

#------------------------------------------------------------
#   Water Tests
#------------------------------------------------------------

def test_water_initialization():
    w = Water(formation_volume_factor=1.0, viscosity=1.0, density=1.0)
    assert w.formation_volume_factor == 1.0
    assert w.viscosity == 1.0
    assert w.density == 1.0

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
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=mock_z_correlation, viscosity_correlation=mock_mu_correlation)
    assert g.critical_pressure == 1000
    assert g.critical_temperature == 500
    assert g.molar_mass == 20
    assert g.z_correlation == mock_z_correlation
    assert g.viscosity_correlation == mock_mu_correlation

def test_gas_z():
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=mock_z_correlation, viscosity_correlation=mock_mu_correlation)
    assert g.z(1000, 500) == 1.0

def test_gas_mu():
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=mock_z_correlation, viscosity_correlation=mock_mu_correlation)
    assert g.viscosity(1000, 500) == 0.01

def test_gas_Bg():
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=mock_z_correlation, viscosity_correlation=mock_mu_correlation)
    assert g.formation_volume_factor(1000, 500) == pytest.approx(0.0283 * 500 / 1000, rel=1e-3)

def test_gas_rho():
    g = Gas(critical_pressure=1000, critical_temperature=500, molar_mass=20, z_correlation=mock_z_correlation, viscosity_correlation=mock_mu_correlation)
    assert g.density(1000, 500) == pytest.approx(144 * 1000 * 20 / 1.0 / 1545.349 / 500, rel=1e-3)