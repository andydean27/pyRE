from pyRE.fluid.correlations import *

class Water:
    """
    Water object class for incompressible water
    
    Parameters
    
    formation_volume_factor : float : Formation volume factor [bbl/STB]
    viscosity : float : Viscosity [cP]
    density : float : Density [g/cm3]
    """

    def __init__(
            self,
            formation_volume_factor: float = 1,
            viscosity: float = 0.75,
            density: float = 1
            ):
        """
        Initialise water object
        
        Arguments
            Bw : float : Formation volume factor [bbl/STB]
            mu : float : Viscosity [cP]
            density : float : Density [g/cm3]
        """
        self.formation_volume_factor = formation_volume_factor
        self.viscosity = viscosity
        self.density = density


class Oil:
    pass

class Gas:
    """
    Gas object class
    
    Parameters
    
    critical_pressure : float : Critical pressure [psia]
    critical_temperature : float : Critical temperature [Rankine]
    molar_mass : float : Molar mass [g/mol]
    z_correlation : function : Z-factor correlation
        options: STANDINGKATZ
    viscosity_correlation : function : Viscosity correlation
        options: LEAGONZALEZEAKIN
    """

    def __init__(
            self,
            critical_pressure: float = 673,
            critical_temperature: float = 343,
            molar_mass: float = 16.04,
            z_correlation = STANDINGKATZ,
            viscosity_correlation = LEEGONZALEZEAKIN
            ):
        """
        Initialise gas object
        
        Arguments
            critical_pressure : float : Critical pressure [psia]
            critical_temperature : float : Critical temperature [Rankine]
            molar_mass : float : Molar mass [g/mol]
            z_correlation (optional) : function : Z-factor correlation
                options: STANDINGKATZ
            viscosity_correlation (optional) : function : Viscosity correlation
                options: LEAGONZALEZEAKIN
        """
        self.critical_pressure = critical_pressure
        self.critical_temperature = critical_temperature
        self.molar_mass = molar_mass
        self.z_correlation = z_correlation
        self.viscosity_correlation = viscosity_correlation

    def z(
            self,
            pressure: float,
            temperature: float,
            ):
        """
        Calculate z-factor

        Arguments
            pressure : float : Pressure [psia]
            temperature : float : Temperature [Rankine]

        Returns
            z : float : z-factor
        """
        return self.z_correlation(pressure, temperature, self)
    
    def viscosity(
            self,
            pressure: float,
            temperature: float,
            ):
        """
        Calculate gas viscosity
        
        Arguments
            pressure : float : Pressure [psia]
            temperature : float : Temperature [Rankine]

        Returns
            viscosity : float : Viscosity [cP]
        """
        
        return self.viscosity_correlation(pressure, temperature, self)
    
    def formation_volume_factor(
            self,
            pressure: float,
            temperature: float,
            ):
        """
        Calculate gas formation volume factor

        Arguments
            pressure : float : Pressure [psia]
            temperature : float : Temperature [Rankine]

        Returns
            Bg : float : Formation volume factor [bbl/SCF]
        """
        
        return 0.0283*self.z(pressure, temperature)*temperature/pressure
    
    def density(
            self,
            pressure: float,
            temperature: float,
            ):
        """
        Calculate gas density
        
        Arguments
            pressure : float : Pressure [psia]
            temperature : float : Temperature [Rankine]

        Returns
            density : float : Density [g/cm3]
        """
        return 144 * pressure * self.molar_mass / self.z(pressure, temperature) / 1545.349 / temperature
    
    def compressibility(
            self,
            pressure: float,
            temperature: float,
            ):
        """
        Calculate gas compressibility
        
        Arguments
            P : float : Pressure [psia]
            T : float : Temperature [Rankine]

        Returns
            Compressibility : float : Compressibility [1/psi]
              = 1/p - 1/z (dz/dp)_T
              = 1/p for pressures under 3000psia
        """
        if pressure <= 3000:
            return 1/pressure

        dp = 0.01
        derivative = (self.z(pressure + 0.01, temperature) - self.z(pressure - 0.01, temperature)) / (2*dp)
        
        return 1/pressure - 1/self.z(pressure, temperature) * derivative
    

    def pseudo_pressure(
            self,
            pressure: float,
            temperature: float,
            ):
        """
        TODO:
        Calculate pseudopressure
        
        Arguments
            P : float : Pressure [psia]
            T : float : Temperature [Rankine]

        Returns
            Pseudo-pressure : float : Pseudopressure
        """
        pass


    
