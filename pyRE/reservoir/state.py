class Saturations:
    """
    Class that holds the saturations of the reservoir
    """
    def __init__(self, **kwargs):
        """
        Initialises the Saturations class
        Arguments:
            **kwargs:
                oil : float : Oil saturation (between 0 and 1)
                gas : float : Gas saturation (between 0 and 1)
                water : float : Water saturation (between 0 and 1)
        """
        self.oil = kwargs.get('oil', None)
        self.gas = kwargs.get('gas', None)
        self.water = kwargs.get('water', None)

        self._validate_saturations()

    def __repr__(self):
        return f"Saturations(oil: {self.oil}, gas: {self.gas}, water: {self.water})"
    
    def _validate_saturations(self):
        """
        Validate that the saturations are within the range [0, 1] and sum to 1.
        """
        if self.oil is not None and (self.oil < 0 or self.oil > 1):
            raise ValueError(f"Oil saturation must be between 0 and 1, got {self.oil}")
        if self.gas is not None and (self.gas < 0 or self.gas > 1):
            raise ValueError(f"Gas saturation must be between 0 and 1, got {self.gas}")
        if self.water is not None and (self.water < 0 or self.water > 1):
            raise ValueError(f"Water saturation must be between 0 and 1, got {self.water}")

        total_saturation = self.sum()
        if total_saturation > 1:
            raise ValueError(f"Total saturation exceeds 1, got {total_saturation}")

    def sum(self):
        """
        Calculate the sum of the saturations.
        """
        return sum(sat for sat in [self.oil, self.gas, self.water] if sat is not None)

class ReservoirState:
    """
    Class that holds the state of the reservoir
    """
    def __init__(self, 
                    pressure: float, 
                    temperature: float,
                    saturations: Saturations):
        """
        Initialises the ReservoirState class
        Arguments:
            pressure : float : Pressure of the reservoir in psia
            temperature : float : Temperature of the reservoir in R
            saturations : Saturations : Saturations object containing oil, gas, and water saturations
        """
        self.pressure = pressure
        self.temperature = temperature
        self.saturations = saturations

    def __repr__(self):
        return f"ReservoirState(pressure: {self.pressure}, temperature: {self.temperature}, saturations: {self.saturations})"
    
    def _validate_state(self):
        """
        Validate the state of the reservoir.
        """
        self.saturations._validate_saturations()

class SimulationState:
    """
    Class that holds the state of the simulation
    """
    def __init__(self, 
                    reservoir_state: ReservoirState,
                    time: float,
                    cumulative_oil_production: float,
                    cumulative_gas_production: float,
                    cumulative_water_production: float):
        """
        Initialises the SimulationState class
        Arguments:
            reservoir_state : ReservoirState : ReservoirState object containing the current state of the reservoir
            time : float : Time in days since the start of the simulation
        """
        self.reservoir_state = reservoir_state
        self.time = time

    def __repr__(self):
        return f"SimulationState(reservoir_state: {self.reservoir_state}, time: {self.time})"
