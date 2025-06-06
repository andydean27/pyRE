
from pyRE.reservoir.models import Reservoir
from pyRE.fluid.models import Oil, Gas, Water
from pyRE.reservoir.state import ReservoirState, Saturations

import math

shape_factors = {
    'circular': 31.62,
    'hexagonal': 31.6,
    'square': 30.8828,
    'equilateral_triangle': 27.6,
}

# Vertical well inflow equations
def _pseudo_steady_flow_rate(
        reservoir: Reservoir,
        state: ReservoirState,
        fluid: str,
        well_flowing_pressure: float,
        wellbore_radius: float,
        skin: float,
        shape_factor: str | float = 'circular',
        pressure_squared_approximation: bool = True) -> float:
    """
    Calculate the pseudo-steady flow rate for a vertical well.
    
    Arguments:
    - reservoir: Reservoir object containing reservoir properties.
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

    if isinstance(shape_factor, str):
        shape_factor = shape_factors.get(shape_factor.lower(), 31.62)
        if shape_factor is None:
            raise ValueError(f"Invalid shape factor: {shape_factor}")
    elif not isinstance(shape_factor, (int, float)):
        raise TypeError(f"Shape factor must be a string or a numeric value, got {type(shape_factor)}")    
    
    # Check if flowing pressure is below reservoir pressure
    if well_flowing_pressure > state.pressure:
        return 0.0


    # Get the fluid model
    if fluid == 'oil':
        if reservoir.fluids.oil is None:
            raise ValueError("Oil model is not defined in the reservoir.")
        fluid = reservoir.fluids.oil
    elif fluid == 'gas':
        if reservoir.fluids.gas is None:
            raise ValueError("Gas model is not defined in the reservoir.")
        fluid = reservoir.fluids.gas
    elif fluid == 'water':
        if reservoir.fluids.water is None:
            raise ValueError("Water model is not defined in the reservoir.")
        fluid = reservoir.fluids.water
    else:
        raise ValueError(f"Invalid fluid type: {fluid}. Must be 'oil', 'gas', or 'water'.")

    # Calculate relative permeabilities
    reservoir.relative_permeability.calculate(state.saturations)
    if fluid == reservoir.fluids.oil:
        relative_permeability = reservoir.relative_permeability.oil
    elif fluid == reservoir.fluids.gas:
        relative_permeability = reservoir.relative_permeability.gas
    elif fluid == reservoir.fluids.water:
        relative_permeability = reservoir.relative_permeability.water
    else:
        raise ValueError(f"Fluid {fluid} not found in the reservoir's relative permeability model.")
    
    # Calculate for gas
    if isinstance(fluid, Gas):
        if pressure_squared_approximation:
            return (
                ((state.pressure**2 - well_flowing_pressure**2) * reservoir.rock.permeability * relative_permeability * reservoir.thickness * reservoir.net_to_gross) /
                (1422 * state.temperature * fluid.z(state.pressure, state.temperature) * fluid.viscosity(state.pressure, state.temperature)) / 
                (0.5 * math.log(10.06 * reservoir.area * 43560 / (shape_factor * wellbore_radius**2)) - 0.75 + skin)
            )
        else:
            raise NotImplementedError("pseudo pressure not implemented yet")
        
    # Calculate for oil or water
    return (
        ((state.pressure - well_flowing_pressure) * reservoir.rock.permeability * relative_permeability * reservoir.thickness * reservoir.net_to_gross) /
        (141.2 * fluid.formation_volume_factor(state.pressure, state.temperature) * fluid.viscosity(state.pressure, state.temperature)) / 
        (0.5 * math.log(10.06 * reservoir.area * 43560 / (shape_factor * wellbore_radius**2)) - 0.75 + skin)
    )