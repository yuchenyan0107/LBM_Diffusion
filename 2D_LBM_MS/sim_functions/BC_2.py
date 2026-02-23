from .common import *
from .boundary_conditions import *


def top_moving_left_intake_bottom_absorb_4_species(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]): # for each direction
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))
    '''
    # Right outlet
    for i in (6, 7, 8):
        g_streamed[:, i, -1, :] = g_dagger_s[:, i, -2, :]
    '''

    # Top wall moving
    '''
    rho_top_wall = 3* lbm_config.c_top[:, None] * xp.ones(lbm_config.nx)[None, :]
    ux, uy = lbm_config.v_top  # wall velocity

    u_lid = xp.array(lbm_config.v_top, dtype=DTYPE)[0]  # 0-D CuPy scalar is OK now
    rho_top_s = 3.0 * lbm_config.c_top[:, None] * xp.ones((1, lbm_config.nx), dtype=DTYPE)

    for i in (5, 4, 6):  # CY=-1 directions (incoming from above)
        c_dot_u = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy  # scalar
        g_streamed[:, i, :, -1] = (
                g_dagger_s[:, OPPOSITE[i], :, -1]
                + 6.0 * w[i] * rho_top_wall * c_dot_u  # since cs^2 = 1/3
        )
    '''
    #g_streamed = top_wall_dirichlet(g_streamed, g_dagger_s, lbm_config)
    g_streamed = top_wall_velocity_dirichlet(g_streamed, lbm_config, x_slice=slice(None))

    # Bottom wall absorption
    g_streamed = bottom_absorption_four_components(g_streamed, g_dagger_s, lbm_config)

    # Left inlet
    #rho_in = 3.0 * lbm_config.c_top[:, None] * xp.ones(lbm_config.ny, dtype=g_dagger_s.dtype)[None, :]  # (ns,ny)

    y = xp.arange(lbm_config.ny, dtype=DTYPE)
    u_wall = xp.array(lbm_config.v_top, dtype=DTYPE)[0]
    ux_in = u_wall * (y / (lbm_config.ny - 1.0))
    uy_in = xp.zeros_like(ux_in)

    g_streamed = zou_he_mixture_velocity_inlet_left(g_streamed, lbm_config, ux_in, uy_in, x0=0, y_slice=slice(None))

    '''
    for s in range(lbm_config.c_top.shape[0]):
        feq_in = eq_single_boundary(rho_in[s], lbm_config.phis[s], ux_in, uy_in)  # (ns,9,ny)
        g_streamed[s, :, 0, :] = feq_in[:, :]
    '''
    # right outlet:
    g_streamed = outlet_zero_gradient_right(g_streamed, lbm_config)


    return g_streamed

def top_wall_dirichlet(g_streamed, g_dagger_s, lbm_config):
    boundary_indexes = xp.array([4,5,6])
    opposite_indexes = OPPOSITE[boundary_indexes]
    row_index = -1

    for s in range(len(lbm_config.c_top)):
        C_w = xp.ones(lbm_config.nx, dtype=DTYPE) * lbm_config.c_top[s] * 3
        # velocity for f_eq here must be 0 for this boundary condition
        #ux_array = lbm_config.v_top[0] * xp.ones(lbm_config.nx, dtype=DTYPE) # doesn't work at all
        ux_array = xp.zeros(g_dagger_s.shape[2])
        f_w_bottom = eq_single_boundary(C_w, lbm_config.phis[s], ux_array,
                                        xp.zeros(g_dagger_s.shape[2]))
        # apply boundary concentration C_w to g_streamed
        g_streamed[s, boundary_indexes, :, row_index] = -g_streamed[s, opposite_indexes, :, row_index] + 2 * f_w_bottom[boundary_indexes,:]

    return g_streamed

def top_wall_velocity_dirichlet(g_streamed, lbm_config, x_slice=slice(None)):
    incoming = xp.array([4, 5, 6])              # your CY=-1 set
    opp = OPPOSITE[incoming]
    row_index = -1
    nx = lbm_config.nx

    ux = xp.full(nx, lbm_config.v_top[0], dtype=g_streamed.dtype)
    uy = xp.full(nx, lbm_config.v_top[1], dtype=g_streamed.dtype)

    # freeze boundary plane before overwriting
    g_pre = g_streamed[:, :, :, row_index].copy()      # (ns, 9, nx)

    for s in range(len(lbm_config.c_top)):
        Cw = xp.full(nx, 3.0 * lbm_config.c_top[s], dtype=g_streamed.dtype)

        feq = eq_single_boundary(Cw, lbm_config.phis[s], ux, uy)  # (9, nx)

        # NEE (bounce-back of non-equilibrium)
        g_streamed[s, incoming, x_slice, row_index] = (
             + feq[incoming, x_slice] + g_streamed[s, opp, :, row_index] - feq[opp, x_slice]
        )

        # enforce exact Dirichlet density by adjusting rest population
        C_now = xp.sum(g_streamed[s, :, :, row_index], axis=0)           # (nx,)
        g_streamed[s, 0, x_slice, row_index] += (Cw[x_slice] - C_now[x_slice])

    return g_streamed

def bottom_absorption_four_components(g_streamed, g_dagger_s, lbm_config):
    boundary_indexes = xp.array([1, 2, 8])
    opposite_indexes = OPPOSITE[boundary_indexes]
    mask_x = (lbm_config.non_absorb_mask == 1)
    mask_idx = xp.nonzero(mask_x)[0]
    row_index = 0
    # component 1, H2: reflecting boundary: no-slip wall
    s=0
    g_streamed[s, boundary_indexes, :, row_index] = g_dagger_s[s, opposite_indexes, :, row_index]

    # component 2 TMIn: Robin condition: dC/dx = k C
    s=1
    b1, b2, b3 = lbm_config.comp2_bc
    b3_rho = b3 * 3
    ne_dx = 1
    C_f2 = xp.sum(g_dagger_s[s, :, :, row_index], axis=0)  # last component, bottom row
    C_w = (C_f2 + ne_dx * b3_rho / b1 / 2) / (1 + ne_dx * b2 / b1 / 2)
    # velocity for f_eq here must be 0 for this boundary condition
    f_w_bottom = eq_single_boundary(C_w, lbm_config.phis[1], xp.zeros(g_dagger_s.shape[2]),xp.zeros(g_dagger_s.shape[2]))
    # apply boundary concentration C_w to g_streamed
    g_streamed[s, boundary_indexes, :, row_index] = -g_dagger_s[s, opposite_indexes, :, row_index] + 2 * f_w_bottom[boundary_indexes, :]
    g_streamed[s, boundary_indexes[:, None], mask_idx[None, :], row_index] = g_dagger_s[s, opposite_indexes[:, None], mask_idx[None, :], row_index]

    #grad2 = -(xp.sum(g_dagger_s[s, :, :, row_index+2], axis=0)-xp.sum(g_dagger_s[s, :, :, row_index], axis=0))/2
    #J2 = -lbm_config.D_s[s] * grad2

    J2 = -lbm_config.D_s[s] * (C_w - C_f2) / (ne_dx/2)

    # component 3 PH3: same absorption flux as component 2
    s = 2
    J3 = J2 * lbm_config.absorption_ratio[s]
    C_f3 = xp.sum(g_dagger_s[s, :, :, row_index], axis=0)
    C_w3 = C_f3 - (ne_dx/2) * (J3 / lbm_config.D_s[s])
    feq_w3 = eq_single_boundary(C_w3, lbm_config.phis[s], xp.zeros(g_dagger_s.shape[2]),xp.zeros(g_dagger_s.shape[2]))
    g_streamed[s, boundary_indexes, :, row_index] = -g_dagger_s[s, opposite_indexes, :, row_index] + 2 * feq_w3[boundary_indexes, :]
    g_streamed[s, boundary_indexes[:, None], mask_idx[None, :], row_index] = g_dagger_s[s, opposite_indexes[:, None], mask_idx[None, :], row_index]

    # component 4, CH4: positive & three times slope of PH3 consumption
    s = 3
    J4 = J2 * lbm_config.absorption_ratio[s]
    C_f4 = xp.sum(g_dagger_s[s, :, :, row_index], axis=0)
    C_w4 = C_f4 - (ne_dx / 2) * (J4 / lbm_config.D_s[s])
    feq_w4 = eq_single_boundary(C_w4, lbm_config.phis[s], xp.zeros(g_dagger_s.shape[2]), xp.zeros(g_dagger_s.shape[2]))
    g_streamed[s, boundary_indexes, :, row_index] = -g_dagger_s[s, opposite_indexes, :, row_index] + 2 * feq_w4[boundary_indexes,:]
    g_streamed[s, boundary_indexes[:, None], mask_idx[None, :], row_index] = g_dagger_s[s, opposite_indexes[:, None], mask_idx[None, :], row_index]
    return g_streamed

def top_wall_noslip(g_streamed, g_dagger_s, lbm_config):
    boundary_indexes = xp.array([4,5,6])
    opposite_indexes = OPPOSITE[boundary_indexes]
    row_index = -1
    # component 1, H2: reflecting boundary: no-slip wall
    for s in range(4):
        g_streamed[s, boundary_indexes, :, row_index] = g_dagger_s[s, opposite_indexes, :, row_index]

    return g_streamed


def tube_flow_bottom_absorb(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]): # for each direction
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))
    '''
    # Right outlet
    for i in (6, 7, 8):
        g_streamed[:, i, -1, :] = g_streamed[:, i, -2, :]
    '''
    # Top wall stationary
    rho_top_wall = 3* lbm_config.c_top[:, None] * xp.ones(lbm_config.nx)[None, :]
    ux, uy = 0 ,0
    for i in (5, 4, 6):  # CY=-1 directions (incoming from above)
        c_dot_u = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy  # scalar
        g_streamed[:, i, :, -1] = (
                g_dagger_s[:, OPPOSITE[i], :, -1]
                + 6.0 * w[i] * rho_top_wall * c_dot_u  # since cs^2 = 1/3
        )

    # Bottom wall absorption
    g_streamed = bottom_absorption_four_components(g_streamed, g_dagger_s, lbm_config)

    # Left inlet (Poiseuille / quadratic profile between no-slip walls)
    rho_in = 3.0 * lbm_config.c_top[:, None] * xp.ones(
        lbm_config.ny, dtype=g_dagger_s.dtype
    )[None, :]  # (ns, ny)

    y = xp.arange(lbm_config.ny, dtype=DTYPE)
    H = DTYPE(lbm_config.ny - 1.0)
    eta = y / H  # 0..1

    u_max = xp.array(lbm_config.v_top, dtype=DTYPE)[0]
    ux_in = 4.0 * u_max * eta * (1.0 - eta)

    uy_in = xp.zeros_like(ux_in)

    #C_in_s = 3.0 * lbm_config.c_top[:, None] * xp.ones(lbm_config.ny, dtype=g_dagger_s.dtype)[None, :]

    #g_streamed = inlet_nee_left(g_streamed, g_dagger_s, lbm_config, rho_in, ux_in, uy_in, y_slice=slice(None))
    '''
    for s in range(lbm_config.c_top.shape[0]):
        feq_in = eq_single_boundary(rho_in[s], lbm_config.phis[s], ux_in, uy_in)
        g_streamed[s, :, 0, :] = feq_in[:, :]
    '''
    g_streamed = zou_he_mixture_velocity_inlet_left(g_streamed, lbm_config, ux_in, uy_in, x0=0, y_slice=slice(None))

    # right outlet:
    #g_streamed = outlet_nee_right(g_streamed, g_dagger_s, lbm_config)
    g_streamed = outlet_zero_gradient_right(g_streamed, lbm_config)

    return g_streamed

def zou_he_mixture_velocity_inlet_left(
    g_streamed,         # (ns, 9, nx, ny) after streaming, will be modified in-place
    lbm_config,
    ux_in, uy_in,       # (ny,) desired inlet velocity profile
    x0=0,
    y_slice=slice(1, -1),
    phi_total=1.0,      # use 1.0 -> standard fluid equilibrium
):
    """
    Mixture-level Zou/He velocity inlet at left boundary x=x0.

    Steps:
      1) form total populations f_tot = sum_s g_s at the inlet plane
      2) compute rho_tot(y) via Zou/He using f_tot and prescribed u(y)
      3) split rho_tot into rho_s = x_s * rho_tot using composition fractions from lbm_config.c_top
      4) reconstruct each species' unknown incoming populations:
           g_s[i] = g_s[opp(i)] + feq_s[i](rho_s,u) - feq_s[opp(i)](rho_s,u)
         for directions with CX>0 at the left boundary.
    """

    ns, q, nx, ny = g_streamed.shape
    assert q == 9, "This inlet assumes D2Q9."

    # --- direction sets (sign-based, indexing-agnostic)
    incoming = xp.where(D2Q9_CX > 0)[0]      # unknown at x=0 (coming into domain)
    known    = xp.where(D2Q9_CX <= 0)[0]     # known at x=0 after streaming
    opp      = OPPOSITE[incoming]

    # --- freeze boundary values before we overwrite anything
    # shape: (ns, 9, ny)
    g0_pre = g_streamed[:, :, x0, :].copy()

    # --- mixture total populations at x=0
    # shape: (9, ny)
    f_tot = xp.sum(g0_pre, axis=0)

    # --- build equilibrium "coefficients" B_i(u) s.t. feq_i = rho * B_i
    ones = xp.ones(ny, dtype=g_streamed.dtype)
    B_tot = eq_single_boundary(ones, phi_total, ux_in, uy_in)  # (9, ny)

    # --- Zou/He solve for rho_tot(y)
    # rho = sum_known f + sum_incoming( f_opp + rho*(B_i - B_opp) )
    # => rho*(1 - sum_incoming(B_i - B_opp)) = sum_known f + sum_incoming f_opp
    K = xp.sum(B_tot[incoming, :] - B_tot[opp, :], axis=0)     # (ny,)
    S_known = xp.sum(f_tot[known, :], axis=0)                  # (ny,)
    S_opp   = xp.sum(f_tot[opp,  :], axis=0)                   # (ny,)

    denom = (1.0 - K)
    denom = xp.where(xp.abs(denom) < 1e-12, 1e-12, denom)
    rho_tot = (S_known + S_opp) / denom                        # (ny,)

    # --- split mixture density into species partial densities using c_top fractions
    frac = xp.array(lbm_config.c_top, dtype=g_streamed.dtype)  # (ns,)
    frac_sum = xp.sum(frac)
    frac = frac / (frac_sum if frac_sum != 0 else 1.0)

    # rho_s: (ns, ny)
    rho_s = frac[:, None] * rho_tot[None, :]

    # --- species-wise Zou/He reconstruction for unknown incoming directions
    # Only touch y_slice to avoid fighting top/bottom BCs at corners.
    for s in range(ns):
        feq_s = eq_single_boundary(rho_s[s], lbm_config.phis[s], ux_in, uy_in)  # (9, ny)

        # known opposite populations come from pre-BC boundary values
        g_opp = g0_pre[s, opp, :]  # (n_in, ny)

        # Zou/He reconstruction (non-equilibrium bounce-back form)
        g_in = g_opp + (feq_s[incoming, :] - feq_s[opp, :])  # (n_in, ny)

        g_streamed[s, incoming, x0, y_slice] = g_in[:, y_slice]

    return g_streamed#, rho_tot


def inlet_nee_left(g_streamed, g_dagger_s, lbm_config, C_in_s, ux_in, uy_in, y_slice=slice(None)):
    """
    NEE inlet at x=0.
    - Only sets incoming directions (cx>0) to avoid overwriting outgoing ones.
    - Uses non-equilibrium from interior plane x=1.
    Shapes:
      g_* : (ns, 9, nx, ny)
      C_in_s : (ns, ny)   (species concentration/density at inlet)
      ux_in, uy_in : (ny,)
    """
    x0, x1 = 0, 1
    incoming_left = xp.where(D2Q9_CX > 0)[0]  # directions entering domain from left boundary

    # interior macro (per species) at x=1
    C_1_s = xp.sum(g_dagger_s[:, :, x1, :], axis=1)  # (ns, ny)

    for s in range(g_dagger_s.shape[0]):
        feq_in = eq_single_boundary(C_in_s[s], lbm_config.phis[s], ux_in, uy_in)  # (9, ny)
        feq_1  = eq_single_boundary(C_1_s[s],  lbm_config.phis[s], ux_in, uy_in)  # (9, ny)

        f_1 = g_dagger_s[s, :, x1, :]  # (9, ny)

        g_streamed[s, incoming_left, x0, y_slice] = (
            feq_in[incoming_left, y_slice] + (f_1[incoming_left, y_slice] - feq_1[incoming_left, y_slice])
        )

    return g_streamed

def outlet_zero_gradient_right(g_streamed, lbm_config):
    """
    Simple Zero-Gradient (Open) Outlet at x=nx-1.
    Copies populations from x=nx-2 to x=nx-1 for incoming directions.
    """
    x_out = lbm_config.nx - 1
    x_in = lbm_config.nx - 2

    # Directions entering the domain from the right boundary (moving Left)
    incoming_from_right = xp.where(D2Q9_CX < 0)[0]

    # Simply copy the distribution from the neighbor node
    # This acts as a "conveyor belt" letting structures exit smoothly.
    for s in range(g_streamed.shape[0]):
        g_streamed[s, incoming_from_right, x_out, :] = g_streamed[s, incoming_from_right, x_in, :]

    return g_streamed


def demo2(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]): # for each direction
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    g_streamed = top_wall_velocity_dirichlet(g_streamed, lbm_config, x_slice=slice(None))

    # Bottom wall absorption
    g_streamed = demo2_bottom(g_streamed, g_dagger_s, lbm_config)

    # Left inlet
    #rho_in = 3.0 * lbm_config.c_top[:, None] * xp.ones(lbm_config.ny, dtype=g_dagger_s.dtype)[None, :]  # (ns,ny)

    y = xp.arange(lbm_config.ny, dtype=DTYPE)
    u_wall = xp.array(lbm_config.v_top, dtype=DTYPE)[0]
    ux_in = u_wall * (y / (lbm_config.ny - 1.0))
    uy_in = xp.zeros_like(ux_in)

    g_streamed = zou_he_mixture_velocity_inlet_left(g_streamed, lbm_config, ux_in, uy_in, x0=0, y_slice=slice(None))


    # right outlet:
    g_streamed = outlet_zero_gradient_right(g_streamed, lbm_config)


    return g_streamed

def demo2_bottom(g_streamed, g_dagger_s, lbm_config):
    boundary_indexes = xp.array([1, 2, 8])
    opposite_indexes = OPPOSITE[boundary_indexes]
    mask_x = (lbm_config.non_absorb_mask == 1)
    mask_idx = xp.nonzero(mask_x)[0]
    row_index = 0
    # component 1, H2: reflecting boundary: no-slip wall
    s=0
    g_streamed[s, boundary_indexes, :, row_index] = g_dagger_s[s, opposite_indexes, :, row_index]

    # component 2 TMIn: Robin condition: dC/dx = k C
    s=1
    b1, b2, b3 = lbm_config.comp2_bc
    b3_rho = b3 * 3
    ne_dx = 1
    C_f2 = xp.sum(g_dagger_s[s, :, :, row_index], axis=0)  # last component, bottom row
    C_w = (C_f2 + ne_dx * b3_rho / b1 / 2) / (1 + ne_dx * b2 / b1 / 2)
    # velocity for f_eq here must be 0 for this boundary condition
    f_w_bottom = eq_single_boundary(C_w, lbm_config.phis[1], xp.zeros(g_dagger_s.shape[2]),xp.zeros(g_dagger_s.shape[2]))
    # apply boundary concentration C_w to g_streamed
    g_streamed[s, boundary_indexes, :, row_index] = -g_dagger_s[s, opposite_indexes, :, row_index] + 2 * f_w_bottom[boundary_indexes, :]
    g_streamed[s, boundary_indexes[:, None], mask_idx[None, :], row_index] = g_dagger_s[s, opposite_indexes[:, None], mask_idx[None, :], row_index]

    return g_streamed