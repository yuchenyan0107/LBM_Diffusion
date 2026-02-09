class lbm_parameters:
    def __init__(self,
                 nx,ny,
                 steps,
                 output_interval,
                 molecular_weights,
                 multiplier,
                 phis,
                 D_s = None,
                 non_absorb_mask = None,
                 bc_bottom = None,
                 bc_top = None,
                 v_top = None, # top wall velocity, vx, vy
                 c_top = None, # top wall partial pressure
                 comp2_bc = None,
                 absorption_ratio = None,
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
        self.comp2_bc = comp2_bc
        self.D_s = D_s
        self.absorption_ratio = absorption_ratio