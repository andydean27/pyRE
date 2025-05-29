from scipy.interpolate import interp1d
from pyRE.reservoir.state import Saturations

class RelativePermeability:
    """
    Base class for relative permeability models
    """

    def __init__(self):
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
            phase_names: list,
            phase_minimum_saturations: list,
            phase_exponents: list,
            phase_max_relative_permeabilities: list):
        """
        Initialise the Brooks-Corey relative permeability model.

        Arguments:
        - phase_names (list): List of phase names (e.g., ['oil', 'gas', 'water']).
        - phase_minimum_saturations (list): List of minimum saturations for each phase.
        - phase_exponents (list): List of exponents for each phase.
        - phase_max_relative_permeabilities (list): List of maximum relative permeabilities for each phase.
        """
        
        if len(phase_names) != len(phase_minimum_saturations) or \
              len(phase_names) != len(phase_exponents) or \
                len(phase_names) != len(phase_max_relative_permeabilities):
            raise ValueError("All input lists must have the same length.")

        # Ensure all phase names are unique
        if len(set(phase_names)) != len(phase_names):
            raise ValueError("Phase names must be unique.")

        self.phase_names = phase_names
        self.phase_minimum_saturations = phase_minimum_saturations
        self.phase_exponents = phase_exponents
        self.phase_max_relative_permeabilities = phase_max_relative_permeabilities

        # Set phase indices
        self.phase_indices = {name: i for i, name in enumerate(phase_names)}

        # Initialise relative permeabilities to None
        self.oil = None
        self.gas = None
        self.water = None

    def __repr__(self):
        return (f"{self.__class__.__name__}(phase_names: {self.phase_names}, "
                f"phase_minimum_saturations: {self.phase_minimum_saturations}, "
                f"phase_exponents: {self.phase_exponents}, "
                f"phase_max_relative_permeabilities: {self.phase_max_relative_permeabilities})")
    def add_phase(self, name: str, minimum_saturation: float, exponent: float, max_relative_permeability: float):
        """
        Add a new phase to the Brooks-Corey model.
        
        Arguments:
        - name: str - Name of the phase (e.g., 'oil', 'gas', 'water').
        - minimum_saturation: float - Minimum saturation for the phase.
        - exponent: float - Exponent for the phase.
        - max_relative_permeability: float - Maximum relative permeability for the phase.
        """
        if name in self.phase_indices:
            raise ValueError(f"Phase '{name}' already exists in the model.")
        
        self.phase_names.append(name)
        self.phase_minimum_saturations.append(minimum_saturation)
        self.phase_exponents.append(exponent)
        self.phase_max_relative_permeabilities.append(max_relative_permeability)
        
        # Update phase indices
        self.phase_indices[name] = len(self.phase_names) - 1
    
    def remove_phase(self, name: str):
        """
        Remove a phase from the Brooks-Corey model.
        
        Arguments:
        - name: str - Name of the phase to remove.
        
        Raises:
        - ValueError: If the phase does not exist in the model.
        """
        if name not in self.phase_indices:
            raise ValueError(f"Phase '{name}' does not exist in the model.")
        
        phase_index = self.phase_indices[name]
        
        # Remove phase data
        del self.phases[phase_index]
        del self.phase_minimum_saturations[phase_index]
        del self.phase_exponents[phase_index]
        del self.phase_max_relative_permeabilities[phase_index]
        
        # Update phase indices
        del self.phase_indices[name]
        for i in range(phase_index, len(self.phase_names)):
            self.phase_indices[self.phase_names[i]] = i

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
        
        phase_index = self.phase_indices.get(phase_name)
        if phase_index is None:
            raise ValueError(f"Phase '{phase_name}' not found in the model phases.")

        min_saturation = self.phase_minimum_saturations[phase_index]
        exponent = self.phase_exponents[phase_index]
        max_relative_permeability = self.phase_max_relative_permeabilities[phase_index]

        # Total irreducible saturation for all phases
        total_irreducible_saturation = sum(self.phase_minimum_saturations)

        # Ensure saturation is above minimum
        if saturation < min_saturation:
            return 0.0
        
        # Calculate relative permeability
        return min(1, max_relative_permeability * ((saturation - min_saturation) / (1.0 - total_irreducible_saturation)) ** exponent)

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
                if phase in self.phase_indices:
                    raise ValueError(f"{phase} is not present (saturation None), but relative permeability curves are defined.")
            else:
                if phase not in self.phase_indices:
                    raise ValueError(f"{phase} phase not found in the model phases.")
                phase_index = self.phase_indices[phase]
                setattr(self, phase, self.calculate_phase_relative_permeability(saturation, phase_index))

