class System():
    def __init__(
        self, 
        box : np.ndarray,
        time : int = 0,
        E_pot : float = 0.,
        E_kin : float = 0.
    ):

    self.box = box
    self.time = time
    self.E_pot = E_pot
    self.E_kin = E_kin

    @property
    def E_tot(self):
        return self.E_pot + self.E_kin
        