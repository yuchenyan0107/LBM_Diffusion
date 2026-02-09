from .common import *

def lattice_stream(f, lbm_config):

    f_streamed = xp.zeros_like(f)
    for i in range(w.shape[0]): # for each direction
        f_streamed[:, i, :, :] = xp.roll(f[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    return f_streamed

def top_moving_test(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]): # for each direction
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    rho_top_wall = 3* lbm_config.c_top[:, None] * xp.ones(lbm_config.nx)[None, :]
    ux, uy = lbm_config.v_top  # wall velocity
    for i in (5, 4, 6):  # CY=-1 directions (incoming from above)
        c_dot_u = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy  # scalar
        g_streamed[:, i, :, -1] = (
                g_dagger_s[:, OPPOSITE[i], :, -1]
                + 6.0 * w[i] * rho_top_wall * c_dot_u  # since cs^2 = 1/3
        )

    for i in xp.array([1,2,8]):
        g_streamed[:, i, :, 0] = g_dagger_s[:, OPPOSITE[i], :, 0]

    return g_streamed

def top_moving_left_intake(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]): # for each direction
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    # Right outlet
    for i in (6, 7, 8):
        g_streamed[:, i, -1, 1:-1] = g_streamed[:, i, -2, 1:-1]

    # Left inlet
    rho_in = 3.0 * lbm_config.c_top[:, None] * xp.ones(lbm_config.ny, dtype=g_dagger_s.dtype)[None, :]  # (ns,ny)

    y = xp.arange(lbm_config.ny, dtype=DTYPE)
    u_wall = xp.array(lbm_config.v_top, dtype=DTYPE)[0]
    ux_in = u_wall * (y / (lbm_config.ny - 1.0))
    uy_in = xp.zeros_like(ux_in)

    for s in range(lbm_config.c_top.shape[0]):
        feq_in = eq_single_boundary(rho_in[s], lbm_config.phis[s], ux_in, uy_in)  # (ns,9,ny)
        g_streamed[s, :, 0, 1:-1] = feq_in[:, 1:-1]

    # Top wall moving
    rho_top_wall = 3* lbm_config.c_top[:, None] * xp.ones(lbm_config.nx)[None, :]
    ux, uy = lbm_config.v_top  # wall velocity
    for i in (5, 4, 6):  # CY=-1 directions (incoming from above)
        c_dot_u = D2Q9_CX[i] * ux + D2Q9_CY[i] * uy  # scalar
        g_streamed[:, i, :, -1] = (
                g_dagger_s[:, OPPOSITE[i], :, -1]
                + 6.0 * w[i] * rho_top_wall * c_dot_u  # since cs^2 = 1/3
        )

    # Bottom wall stationary
    for i in xp.array([1,2,8]):
        g_streamed[:, i, :, 0] = g_dagger_s[:, OPPOSITE[i], :, 0]

    return g_streamed

def top_moving_left_intake_bottom_absorb(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]): # for each direction
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    # Right outlet
    for i in (6, 7, 8):
        g_streamed[:, i, -1, 1:-1] = g_streamed[:, i, -2, 1:-1]

    # Left inlet
    rho_in = 3.0 * lbm_config.c_top[:, None] * xp.ones(lbm_config.ny, dtype=g_dagger_s.dtype)[None, :]  # (ns,ny)

    y = xp.arange(lbm_config.ny, dtype=DTYPE)
    u_wall = xp.array(lbm_config.v_top, dtype=DTYPE)[0]
    ux_in = u_wall * (y / (lbm_config.ny - 1.0))
    uy_in = xp.zeros_like(ux_in)

    for s in range(lbm_config.c_top.shape[0]):
        feq_in = eq_single_boundary(rho_in[s], lbm_config.phis[s], ux_in, uy_in)  # (ns,9,ny)
        g_streamed[s, :, 0, 1:-1] = feq_in[:, 1:-1]

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
    b1, b2, b3, reflection_boundary = lbm_config.bc_bottom
    g_streamed = top_bottom_fixed('bottom', g_streamed, g_dagger_s, b1, b2, b3, reflection_boundary, lbm_config)

    return g_streamed



def top_bottom_fixed(direction, g_streamed, g_dagger_s, b1, b2, b3, reflection_boundary, lbm_config):
    '''
    b1 = 0: fixed boundary concentration
    1/b1 = absorption rate
    b2 = 1
    b3 = concentration at stable
    b3: partial pressure input, then times 3 to get rho (number density)
    '''

    b3_rho = 3*b3

    if direction == 'top':
        boundary_indexes = xp.array([4,5,6])
        row_index = -1
    elif direction == 'bottom':
        boundary_indexes = xp.array([1,2,8])
        row_index = 0

    dx = 1
    ne = 1

    for s in range(g_dagger_s.shape[0]):

        if reflection_boundary[s] == 0:  # if boundary is not reflection:

            if b1[s] == 0:  # Dirichlet boundary
                C_w = xp.ones(g_dagger_s.shape[2], dtype=DTYPE) * b3_rho[s] / b2[s]

            else: # absorption
                C_f = xp.sum(g_dagger_s[s, :, :, row_index], axis=0)  # last component, bottom row

                C_w = (C_f + ne * dx * b3_rho[s] / b1[s] / 2) / (1 + ne * dx * b2[s] / b1[s] / 2)

            # equilibrium based on concentration C_w
            f_w_bottom = eq_single_boundary(C_w, lbm_config.phis[s], xp.zeros(g_dagger_s.shape[2]),
                                            xp.zeros(g_dagger_s.shape[2]))

            for i in boundary_indexes: # apply boundary concentration C_w to g_streamed
                g_streamed[s, i, :, row_index] = -g_dagger_s[s, OPPOSITE[i], :, row_index] + 2 * f_w_bottom[i]

            if b1[s] != 0:  # for absorption condition, overwrite masked area to no-slip wall
                for i in boundary_indexes:
                    g_streamed[s, i, (lbm_config.non_absorb_mask == 1), row_index] = g_dagger_s[s, OPPOSITE[i], (lbm_config.non_absorb_mask == 1), row_index]

        else:  # no-slip wall:
            for i in boundary_indexes:
                g_streamed[s, i, :, row_index] = g_dagger_s[s, OPPOSITE[i], :, row_index]

    return g_streamed


def eq_single_boundary(C_w, phi_s, ux_star, uy_star):
    u_sq = ux_star ** 2 + uy_star ** 2  # size of (nx or...y...)
    cu = D2Q9_CX[:, None] * ux_star[None, :] + D2Q9_CY[:, None] * uy_star[None, :]

    f_w_eq = w[:, None] * C_w[None, :] * (
            phi_s
            + 3 * cu
            + 4.5 * cu ** 2
            - 1.5 * u_sq
    )

    f_w_eq[0, :] = w[0] * C_w * ((9 - 5 * phi_s) / 4 - 1.5 * u_sq)

    return f_w_eq # (9,nx)

def absorb_BC_fix_top_buttom(g_dagger_s, lbm_config):

    g_streamed = xp.zeros_like(g_dagger_s)
    for i in range(w.shape[0]):
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    #bottom absorption:
    b1, b2, b3, reflection_boundary = lbm_config.bc_bottom

    g_streamed = top_bottom_fixed('bottom', g_streamed, g_dagger_s, b1, b2, b3, reflection_boundary, lbm_config)

    # top:
    b1, b2, b3, reflection_boundary = lbm_config.bc_top

    g_streamed = top_bottom_fixed('top', g_streamed, g_dagger_s, b1, b2, b3, reflection_boundary, lbm_config)

    return g_streamed
