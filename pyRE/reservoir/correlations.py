


class Corey2PhaseRelativePermeability:
    def __init__(
            self,
            s1cr: float = 0.35,
            s2cr: float = 0,
            n1: float = 3,
            n2: float = 1,
            kr12c: float = 1,
            kr21c: float = 0.7
    ):
        """
        Initialize the Corey 2-phase relative permeability model.

        Arguments:
        - s1cr (float): Critical saturation of phase 1.
        - s2cr (float): Critical saturation of phase 2.
        - n1 (float): Corey exponent for phase 1.
        - n2 (float): Corey exponent for phase 2.
        - kr12c (float): Relative permeability of phase 1 at critical saturation of phase 2.
        - kr21c (float): Relative permeability of phase 2 at critical saturation of phase 1.
        """
        self.s1cr = s1cr
        self.s2cr = s2cr
        self.n1 = n1
        self.n2 = n2
        self.kr12c = kr12c
        self.kr21c = kr21c
    
    def function(self, s1: float):

        # Validate s1 in range
        s1 = max(self.s1cr, min(s1, 1.0 - self.s2cr))

        # Calculate kr1 and kr2
        kr1 = max(0, min(self.kr12c, self.kr12c * ((s1 - self.s1cr)/(1 - self.s1cr)) ** self.n1))
        kr2 = max(0, min(self.kr21c, self.kr21c * ((1 - s1 - self.s2cr)/(1 - self.s1cr - self.s2cr)) ** self.n2))
        
        return kr1, kr2
