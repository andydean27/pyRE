import math

# Z-factor correlations

def STANDINGKATZ(P: float, T: float, gas_obj):
    """
    Standing-Katz Z-factor correlation
    
    Arguments
        P : float : Pressure [psia]
        T : float : Temperature [Rankine]
        gas_obj : gas : Gas model
    """

    # Calculate pseudoreduced pressure and temperature
    Pr = P / gas_obj.Pc
    Tr = T / gas_obj.Tc

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

def LEEGONZALEZEAKIN(P: float, T: float, gas_obj):
    """
    Lee-Gonzalez-Eakin gas viscosity correlation
    
    Arguments
        P : float : Pressure [psia]
        T : float : Temperature [Rankine]
        gas_obj : gas : Gas model
    """
    x = 2.57 + 1914.5 / T + 0.0095 * gas_obj.M
    Y = 1.11 + 0.04 * x
    k = (7.77 + 0.0063 * gas_obj.M) * T ** 1.5 / (122.4 + 12.8 * gas_obj.M + T)


    # Calculate gas viscosity
    mu = k * math.exp(x * ((gas_obj.rho(P, T) * 0.016) ** Y)) * 0.0001

    return mu