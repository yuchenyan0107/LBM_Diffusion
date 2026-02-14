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
    rho_top_wall = 3* lbm_config.c_top[:, None] * xp.ones(lbm_config.nx)[None, :]
    ux, uy = lbm_config.v_top  # wall velocity
    for i in (5, 4, 6):  # CY=-1 directions (incoming from above)
        c_dot_u = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy  # scalar
        g_streamed[:, i, :, -1] = (
                g_dagger_s[:, OPPOSITE[i], :, -1]
                + 6.0 * w[i] * rho_top_wall * c_dot_u  # since cs^2 = 1/3
        )

    # Bottom wall absorption
    g_streamed = bottom_absorption_four_components(g_streamed, g_dagger_s, lbm_config)

    # Left inlet
    rho_in = 3.0 * lbm_config.c_top[:, None] * xp.ones(lbm_config.ny, dtype=g_dagger_s.dtype)[None, :]  # (ns,ny)

    y = xp.arange(lbm_config.ny, dtype=DTYPE)
    u_wall = xp.array(lbm_config.v_top, dtype=DTYPE)[0]
    ux_in = u_wall * (y / (lbm_config.ny - 1.0))
    uy_in = xp.zeros_like(ux_in)

    g_streamed = inlet_nee_left(g_streamed, g_dagger_s, lbm_config, rho_in, ux_in, uy_in, y_slice=slice(1, -1))

    '''
    for s in range(lbm_config.c_top.shape[0]):
        feq_in = eq_single_boundary(rho_in[s], lbm_config.phis[s], ux_in, uy_in)  # (ns,9,ny)
        g_streamed[s, :, 0, :] = feq_in[:, :]
    '''
    # right outlet:
    g_streamed = outlet_nee_right(g_streamed, g_dagger_s, lbm_config)


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

    g_streamed = inlet_nee_left(g_streamed, g_dagger_s, lbm_config, rho_in, ux_in, uy_in, y_slice=slice(1, -1))
    '''
    for s in range(lbm_config.c_top.shape[0]):
        feq_in = eq_single_boundary(rho_in[s], lbm_config.phis[s], ux_in, uy_in)
        g_streamed[s, :, 0, :] = feq_in[:, :]
    '''
    # right outlet:
    g_streamed = outlet_nee_right(g_streamed, g_dagger_s, lbm_config)

    return g_streamed

def inlet_nee_left(g_streamed, g_dagger_s, lbm_config, C_in_s, ux_in, uy_in, y_slice=slice(1, -1)):
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


def outlet_nee_right(g_streamed, g_dagger_s, lbm_config):
    x_out = lbm_config.nx - 1
    x_in  = lbm_config.nx - 2
    incoming = xp.where(D2Q9_CX < 0)[0]

    # macros at interior plane
    C_in = xp.sum(g_dagger_s[:, :, x_in, :], axis=1)   # (ns, ny)
    # choose outlet macros (simple extrapolation)
    C_out = C_in.copy()

    # if you have mixture velocity field available, use it here
    # otherwise for species-only advection you may use a prescribed u_out or extrapolate from interior
    ux_out = xp.zeros(lbm_config.ny, dtype=g_dagger_s.dtype)
    uy_out = xp.zeros_like(ux_out)

    for s in range(g_dagger_s.shape[0]):
        feq_in  = eq_single_boundary(C_in[s],  lbm_config.phis[s], ux_out, uy_out)   # (9, ny)
        feq_out = eq_single_boundary(C_out[s], lbm_config.phis[s], ux_out, uy_out)

        f_in  = g_dagger_s[s, :, x_in, :]   # (9, ny)
        f_out = feq_out + (f_in - feq_in)

        g_streamed[s, incoming, x_out, 1:-1] = f_out[incoming, 1:-1]

    return g_streamed

