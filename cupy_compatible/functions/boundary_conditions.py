from .use_cupy import *

D2Q9_CX = xp.array([0, 0, 1, 1, 1, 0, -1, -1, -1], dtype=xp.float32) # CW: up, ur, r, dr, d.....
D2Q9_CY = xp.array([0, 1, 1, 0, -1, -1, -1, 0, 1], dtype=xp.float32)

w = xp.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36], dtype=xp.float32)
OPPOSITE = xp.array([0,5,6,7,8,1,2,3,4], dtype=xp.float32)
theta = 0.5

def lattice_stream(f, phi, step, non_absorb_mask, bc_top, bc_bottom):

    f_streamed = xp.zeros_like(f, dtype=xp.float32)
    for i in range(w.shape[0]): # for each species
        f_streamed[:, i, :, :] = xp.roll(f[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    return f_streamed

def lattice_stream_object_periodical(f, phi, step, non_absorb_mask, bc_top, bc_bottom):

    f_streamed = xp.zeros_like(f, dtype=xp.float32)
    for i in range(w.shape[0]): # for each species
        f_streamed[:, i, :, :] = xp.roll(f[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    nx = f_streamed.shape[2]
    ny = f_streamed.shape[3]
    X, Y = xp.meshgrid(xp.arange(nx), xp.arange(ny))
    cylinder = ((X - (nx / 5)) ** 2 + (Y - (ny / 2)) ** 2) < (ny / 8) ** 2
    cylinder = cylinder.T

    inside_boundary = f_streamed[:, :, cylinder]
    inside_boundary = inside_boundary[:, [0, 5, 6, 7, 8, 1, 2, 3, 4], :]  # OPPOSITE
    f_streamed[:, :, cylinder] = inside_boundary


    return f_streamed

def lattice_stream_BC_full(g_dagger_s, phi, step, non_absorb_mask, bc_top, bc_bottom):

    g_streamed = xp.zeros_like(g_dagger_s, dtype=xp.float32)
    for i in range(w.shape[0]):
        g_streamed[:, i, :, :] = xp.roll(g_dagger_s[:, i, :, :], (int(D2Q9_CX[i]), int(D2Q9_CY[i])), axis=(1, 2))

    #bottom absorption:
    b1, b2, b3, reflection_boundary = bc_bottom

    g_streamed = top_bottom_boundary('bottom', g_streamed, g_dagger_s, phi, b1, b2, b3, reflection_boundary, non_absorb_mask)

    # top:
    b1, b2, b3, reflection_boundary = bc_top

    g_streamed = top_bottom_boundary('top', g_streamed, g_dagger_s, phi, b1, b2, b3, reflection_boundary, non_absorb_mask)

    '''
    if step % 200 == 0:
        molecular_weight = xp.array([28.0, 2.0, 44.0])
        nB = 2
        rho_s, ux_s, uy_s, rho_mix, p_mix = calculate_moment(g_streamed, phi)
        m_mix = calculate_m_mix(rho_s, rho_mix, molecular_weight)

        CHI_sc = calculate_CHI(m_mix, molecular_weight, nB)

        ux_star_s, uy_star_s = calculate_u_star(CHI_sc, rho_s, rho_mix, ux_s, uy_s)

        #plot_vector(cp.asnumpy(ux_star_s[0, :, :]), cp.asnumpy(uy_star_s[0, :, :]), zoom = 2)
        #plt.plot(cp.asnumpy(ux_star_s[0, :, :])[:, 200])
    '''
    return g_streamed

def top_bottom_boundary(direction, g_streamed, g_dagger_s, phi, b1, b2, b3, reflection_boundary, non_absorb_mask):

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
                C_w = xp.ones(g_dagger_s.shape[2], dtype=xp.float32) * b3[s] / b2[s]

            else: # absorption
                C_f = xp.sum(g_dagger_s[s, :, :, row_index], axis=0, dtype=xp.float32)  # last component, bottom row

                C_w = (C_f + ne * dx * b3[s] / b1[s] / 2) / (1 + ne * dx * b2[s] / b1[s] / 2)

            # equilibrium based on concentration
            f_w_bottom = eq_single_boundary(C_w, phi[s], xp.zeros(g_dagger_s.shape[2]),
                                            xp.zeros(g_dagger_s.shape[2]))

            for i in boundary_indexes:
                g_streamed[s, i, :, row_index] = -g_dagger_s[s, OPPOSITE[i], :, row_index] + 2 * f_w_bottom[i]

            if b1[s] != 0:  # for absorption condition, overwrite masked area to no-slip wall
                #nx, a = g_dagger_s.shape[2], g_dagger_s.shape[2] // 5
                #non_absorbing = (xp.arange(nx) // a) % 2

                for i in boundary_indexes:
                    g_streamed[s, i, (non_absorb_mask == 1), row_index] = g_dagger_s[s, OPPOSITE[i], (non_absorb_mask == 1), row_index]

        else:  # no-slip wall:
            for i in boundary_indexes:
                g_streamed[s, i, :, row_index] = g_dagger_s[s, OPPOSITE[i], :, row_index]

    return g_streamed


def eq_single_boundary(C_w, phi_s, ux_star, uy_star):
    u_sq = ux_star ** 2 + uy_star ** 2  # (nx or...y...)
    cu = D2Q9_CX[:, None] * ux_star[None, :] + D2Q9_CY[:, None] * uy_star[None, :]

    f_w_eq = w[:, None] * C_w[None, :] * (
            phi_s
            + 3 * cu
            + 4.5 * cu ** 2
            - 1.5 * u_sq
    )

    f_w_eq[0, :] = w[0] * C_w * ((9 - 5 * phi_s) / 4 - 1.5 * u_sq)

    return f_w_eq
