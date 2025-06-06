

class Rock:
    """
    Base class for rock properties.
    """
    def __init__(
            self, 
            porosity: float, 
            permeability: float, 
            density: float,
            compressibility: float
            ):
        """
        Initialises the Rock class with basic properties.
        Arguments:
            - porosity : float : Porosity of the rock (fraction)
            - permeability : float : Permeability of the rock (mD)
            - density : float : Density of the rock (g/cm3)
            - compressibility : float : Compressibility of the rock (1/psi)
        """
        self.porosity = porosity
        self.permeability = permeability
        self.density = density
        self.compressibility = compressibility


class Conventional(Rock):
    """
    Rock properties for conventional reservoirs.
    """
    def __init__(
            self, 
            *args,
            **kwargs
            ):
        super().__init__(*args, **kwargs)


class Coal(Rock):
    """
    Rock properties for coal seam reservoirs.
    This class inherits from the Rock class and adds specific properties for coal.

    """
    def __init__(
            self, 
            gas_content: float,
            langmuir_pressure: float,
            langmuir_volume: float,
            ash_fraction: float = 0,
            moisture_fraction: float = 0,
            *args,
            **kwargs, 
            ):
        """
        Initialises the Coal class with specific properties for coal seam reservoirs.

        Arguments
            - porosity : float : Porosity of the rock (fraction)
            - permeability : float : Permeability of the rock (mD)
            - density : float : Density of the rock (g/cm3)
            - compressibility : float : Compressibility of the rock (1/psi)
            - gas_content : float : Gas content in scf/ton
            - langmuir_pressure : float : Langmuir pressure in psia
            - langmuir_volume : float : Langmuir volume in scf/ton (dry ash-free basis)
            - ash_fraction : float : Ash fraction (default is 0)
            - moisture_fraction : float : Moisture fraction (default is 0)
        """
        super().__init__(*args, **kwargs)

        self.gas_content = gas_content  # Gas content in scf/ton
        self.langmuir_pressure = langmuir_pressure
        self.langmuir_volume = langmuir_volume
        self.ash_fraction = ash_fraction
        self.moisture_fraction = moisture_fraction

        # Calculate critical desorption pressure
        self.critical_desorption_pressure = self.langmuir_pressure/(self.langmuir_volume/self.gas_content - 1)