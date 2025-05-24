from pyRE.fluid.models import Oil, Gas, Water

class FluidCollection:
    """
    Collection of fluid models for the use in reservoir models
    """

    def __init__(self, **kwargs):
        """
        Initialises the FluidCollection

        Arguments:
            **kwargs:
                oil : Oil : Oil model
                gas : Gas : Gas model
                water : Water : Water model
        """

        self.oil = None
        self.gas = None
        self.water = None

        oil = kwargs.get('oil', None)
        gas = kwargs.get('gas', None)
        water = kwargs.get('water', None)

        #check if model is set and is the correct type
        if oil is not None:
            if isinstance(oil, Oil):
                self.oil = oil
            else:
                raise TypeError(f"Expected Oil model, got {type(oil)}")
            
        if gas is not None:
            if isinstance(gas, Gas):
                self.gas = gas
            else:
                raise TypeError(f"Expected Gas model, got {type(gas)}")
            
        if water is not None:
            if isinstance(water, Water):
                self.water = water
            else:
                raise TypeError(f"Expected Water model, got {type(water)}")
            
    def __repr__(self):
        return f"FluidCollection(oil: {self.oil}, gas: {self.gas}, water: {self.water})"
    
    def add_fluid(self, fluid):
        """
        Add a fluid model to the collection.
        """
        if isinstance(fluid, Oil):
            self.oil = fluid
        elif isinstance(fluid, Gas):
            self.gas = fluid
        elif isinstance(fluid, Water):
            self.water = fluid
        else:
            raise TypeError(f"Expected Oil, Gas or Water model, got {type(fluid)}")