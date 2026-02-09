class lbm_parameters:
    def __init__(self,
                 nx,ny,
                 steps,
                 output_interval,
                 molecular_weights,
                 multiplier,
                 phis,
                 non_absorb_mask = None,
                 bc_bottom = None,
                 bc_top = None,
                 v_top = None, # top wall velocity, vx, vy
                 c_top = None, # top wall partial pressure
                 ):
        self.nx = nx
        self.ny = ny
        self.steps = steps
        self.output_interval = output_interval
        self.molecular_weights = molecular_weights
        self.multiplier = multiplier
        self.phis = phis
        self.non_absorb_mask = non_absorb_mask
        self.bc_bottom = bc_bottom
        self.bc_top = bc_top
        self.v_top = v_top
        self.c_top = c_top