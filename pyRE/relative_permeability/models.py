from scipy.interpolate import interp1d
from pyRE.reservoir.state import Saturations
from dataclasses import dataclass

class RelativePermeability:
    """
    Base class for relative permeability models
    """

    def __init__(
            self,
            oil_curve: dict = None,
            gas_curve: dict = None,
            water_curve: dict = None):
        """
        Default behaviour returns a relative permeability of 1 for all phases
        """

    def __repr__(self):
        return f"{self.__class__.__name__}()"
    
    def calculate(self, saturations: Saturations):
        """
        Calculate relative permeability for given saturations.
        
        This method should be overridden by subclasses.
        
        Arguments:
        - saturations: Saturations object containing oil, gas, and water saturations.
        
        Returns:
        - None
        """
        raise NotImplementedError("This method should be overridden by subclasses.")
    
class BrooksCoreyRelativePermeability(RelativePermeability):
    """
    Brooks-Corey relative permeability model.

    Attributes:
    - oil (float): Calculated relative permeability for the oil phase.
    - gas (float): Calculated relative permeability for the gas phase.
    - water (float): Calculated relative permeability for the water phase.
    """

    def __init__(
            self,
            oil_curve: 'BrooksCoreyRelativePermeability.Curve' = None,
            gas_curve: 'BrooksCoreyRelativePermeability.Curve' = None,
            water_curve: 'BrooksCoreyRelativePermeability.Curve' = None):
        """
        Initialise the Brooks-Corey relative permeability model.

        Arguments:
        - oil_curve: Curve object for the oil phase.
        - gas_curve: Curve object for the gas phase.
        - water_curve: Curve object for the water phase.
        Each curve should contain:
        - minimum_saturation (float): Minimum saturation for the phase.
        - exponent (float): Exponent for the relative permeability calculation.
        - max_relative_permeability (float): Maximum relative permeability for the phase.
        """

        self.oil_curve = oil_curve or None
        self.gas_curve = gas_curve or None
        self.water_curve = water_curve or None

        # Calculate total irreducible saturation
        self.total_irreducible_saturation = sum(
            curve.minimum_saturation for curve in [self.oil_curve, self.gas_curve, self.water_curve] if curve is not None
        )

        # Initialise relative permeabilities to None
        self.oil = None
        self.gas = None
        self.water = None

    def __repr__(self):
        return f"{self.__class__.__name__}(oil_curve: {self.oil_curve}, gas_curve: {self.gas_curve}, water_curve: {self.water_curve})"
    
    @dataclass
    class Curve:
        """
        Data class to hold the parameters for a relative permeability curve.
        
        Attributes:
        - minimum_saturation (float): Minimum saturation for the phase.
        - exponent (float): Exponent for the relative permeability calculation.
        - max_relative_permeability (float): Maximum relative permeability for the phase.
        """
        minimum_saturation: float
        exponent: float
        max_relative_permeability: float
        def __repr__(self):
            return f"Curve(minimum_saturation: {self.minimum_saturation}, exponent: {self.exponent}, max_relative_permeability: {self.max_relative_permeability})"

    def calculate_phase_relative_permeability(self, phase_name: str, saturation: float):
        """
        Calculate relative permeability for a single phase based on its saturation.
        
        Arguments:
        - phase_name: str - Name of the phase (e.g., 'oil', 'gas', 'water').
        - saturation: float - Saturation of the phase.

        
        Returns:
        - float - Relative permeability for the phase. Returns 0.0 if saturation is below the minimum.
        """
        # Check if saturation is not None
        if saturation is None:
            return None
        
        # Get the phase curve parameters
        curve = getattr(self, f"{phase_name}_curve", None)
        if curve is None:
            raise ValueError(f"No curve defined for phase '{phase_name}'.")

        min_saturation = curve.minimum_saturation
        exponent = curve.exponent
        max_relative_permeability = curve.max_relative_permeability

        # Ensure saturation is above minimum
        if saturation < min_saturation:
            return 0.0
        
        # Calculate relative permeability
        return min(1, max_relative_permeability * ((saturation - min_saturation) / (1.0 - self.total_irreducible_saturation)) ** exponent)

    def calculate(self, saturations: Saturations):
        """
        Calculate relative permeability for given saturations.
        
        Arguments:
        - saturations: Saturations object containing oil, gas, and water saturations.

        Returns:
        - None
        """

        # Reset relative permeabilities
        self.oil = None
        self.gas = None
        self.water = None

        for phase, saturation in zip(['oil', 'gas', 'water'], [saturations.oil, saturations.gas, saturations.water]):
            if saturation is None:
                # Check if a curve is defined for the phase
                if getattr(self, f"{phase}_curve", None) is not None:
                    raise ValueError(f"{phase} is not present (saturation None), but relative permeability curves are defined.")
            else:

                setattr(self, phase, self.calculate_phase_relative_permeability(phase, saturation))

