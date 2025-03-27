from typing import Optional

class VerticalWell:
    """
    """

    def __init__ (
            self,
            wellbore_radius: float,
            measured_depth: float,
            surface_elevation: float = 0
    ):
        self.wellbore_radius = wellbore_radius
        self.measured_depth = measured_depth
        self.surface_elevation = surface_elevation