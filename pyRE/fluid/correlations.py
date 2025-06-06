import math


# Z-factor correlations

def StandingKatz(P: float, T: float, **kwargs):
    """
    Standing-Katz Z-factor correlation
    
    Arguments
        - P : float : Pressure [psia]
        - T : float : Temperature [Rankine]
        - **kwargs : Keyword arguments
            - critical_pressure : float : Critical pressure [psia]
            - critical_temperature : float : Critical temperature [Rankine]

    Returns:
        - Z : float : Z-factor
    """

    # Keyword arguments
    critical_pressure = kwargs.get('critical_pressure')
    critical_temperature = kwargs.get('critical_temperature')
    if critical_pressure is None or critical_temperature is None:
        raise ValueError("Both critical_pressure and critical_temperature must be provided for Z-factor calculation.")

    # Calculate pseudoreduced pressure and temperature
    Pr = P / critical_pressure
    Tr = T / critical_temperature

    # Standing-Katz correlation
    a = 1.39*(Tr-0.92)**0.5 - 0.36*Tr - 0.101
    e = 9*(Tr-1)
    b = (0.62 - 0.23*Tr)*Pr + (0.066/(Tr-0.86) - 0.037)*Pr**2 + 0.32*Pr**6/(10**e)
    c = 0.132 - 0.32*math.log(Tr)/math.log(10)
    f = 0.3106 - 0.49*Tr + 0.1824*Tr**2
    d = 10**f

    Z = a + (1 - a)*math.exp(-b) + c*Pr**d

    return Z

def DRANCHUKABOUKHAIF(P: float, T: float, Pc: float, Tc: float):
    pass

# Viscosity

def LeeGonzalezEakin(P: float, T: float, **kwargs):
    """
    Lee-Gonzalez-Eakin gas viscosity correlation
    
    Arguments:
        - P : float : Pressure [psia]
        - T : float : Temperature [Rankine]
        - **kwargs : Keyword arguments
            - molar_mass : float : Molar mass [g/mol]
            - density : float : Density 
    """

    # Keyword arguments
    molar_mass = kwargs.get('molar_mass')
    density = kwargs.get('density')
    if molar_mass is None or density is None:
        raise ValueError("Both molar_mass and density must be provided for viscosity calculation.")

    x = 2.57 + 1914.5 / T + 0.0095 * molar_mass
    Y = 1.11 + 0.04 * x
    k = (7.77 + 0.0063 * molar_mass) * T ** 1.5 / (122.4 + 12.8 * molar_mass + T)


    # Calculate gas viscosity
    mu = k * math.exp(x * ((density * 0.016) ** Y)) * 0.0001

    return mu