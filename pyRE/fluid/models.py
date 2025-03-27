from pyRE.fluid.correlations import *
import matplotlib.pyplot as plt

class Water:
    """
    Water object class for incompressible water
    
    Parameters
    
    Bw : float : Formation volume factor [bbl/STB]
    mu : float : Viscosity [cP]
    rho : float : Density [g/cm3]
    """

    def __init__(
            self,
            Bw: float = 1,
            mu: float = 0.75,
            rho: float = 1
            ):
        """
        Initialise water object
        
        Arguments
            Bw : float : Formation volume factor [bbl/STB]
            mu : float : Viscosity [cP]
            rho : float : Density [g/cm3]
        """
        self.Bw = Bw
        self.mu = mu
        self.rho = rho
    
    def load_from_json():
        pass

class Oil:
    pass

class Gas:
    """
    Gas object class
    
    Parameters
    
    Pc : float : Critical pressure [psia]
    Tc : float : Critical temperature [Rankine]
    M : float : Molar mass [g/mol]
    z_correlation : function : Z-factor correlation
        options: STANDINGKATZ
    mu_correlation : function : Viscosity correlation
        options: LEAGONZALEZEAKIN
    """

    def __init__(
            self,
            Pc: float = 673,
            Tc: float = 343,
            M: float = 16.04,
            z_correlation = STANDINGKATZ,
            mu_correlation = LEEGONZALEZEAKIN
            ):
        """
        Initialise gas object
        
        Arguments
            Pc : float : Critical pressure [psia]
            Tc : float : Critical temperature [Rankine]
            M : float : Molar mass [g/mol]
            z_correlation (optional) : function : Z-factor correlation
                options: STANDINGKATZ
            mu_correlation (optional) : function : Viscosity correlation
                options: LEAGONZALEZEAKIN
        """
        self.Pc = Pc
        self.Tc = Tc
        self.M = M
        self.z_correlation = z_correlation
        self.mu_correlation = mu_correlation

    def z(
            self,
            P: float,
            T: float,
            ):
        """
        Calculate z-factor

        Arguments
            P : float : Pressure [psia]
            T : float : Temperature [Rankine]

        Returns
            z : float : z-factor
        """
        return self.z_correlation(P, T, self)
    
    def mu(
            self,
            P: float,
            T: float,
            ):
        """
        Calculate gas viscosity
        
        Arguments
            P : float : Pressure [psia]
            T : float : Temperature [Rankine]

        Returns
            mu : float : Viscosity [cP]
        """
        
        return self.mu_correlation(P, T, self)
    
    def Bg(
            self,
            P: float,
            T: float,
            ):
        """
        Calculate gas formation volume factor

        Arguments
            P : float : Pressure [psia]
            T : float : Temperature [Rankine]

        Returns
            Bg : float : Formation volume factor [bbl/SCF]
        """
        
        return 0.0283*self.z(P, T)*T/P
    
    def rho(
            self,
            P: float,
            T: float,
            ):
        """
        Calculate gas density
        
        Arguments
            P : float : Pressure [psia]
            T : float : Temperature [Rankine]

        Returns
            rho : float : Density [g/cm3]
        """
        return 144 * P * self.M / self.z(P, T) / 1545.349 / T
    
    def load_from_json():
        pass
    
