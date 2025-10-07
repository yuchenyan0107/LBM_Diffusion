import cupy as cp
import numpy as np

def flow_drift_CUDA(F, component_x, component_y):

    '''F = cp.asarray(F)
    component_x = cp.asarray(component_x)
    component_y = cp.asarray(component_y)'''

    N, nx, ny = F.shape

    idx_x = cp.arange(nx)
    idx_y = cp.arange(ny)

    shifted_x = cp.array(idx_x[None, :, None] - component_x[:, None, None]) % nx
    shifted_y = cp.array(idx_y[None, None, :] - component_y[:, None, None]) % ny

    F = F[cp.arange(N)[:, None, None], shifted_x, shifted_y]

    #F = cp.asnumpy(F)

    return F

def inlet_x_CUDA(F, component_x, component_y, w, ux_in, uy_in, rho0, ny):

    '''component_x = cp.asarray(component_x)
    component_y = cp.asarray(component_y)
    w = cp.asarray(w)
    F = cp.asarray(F)'''

    ux_in = cp.full((F.shape[0], ny), ux_in)
    uy_in = cp.full((F.shape[0], ny), uy_in)

    F_inlet_eq = w[:, None] * rho0 * (
        1
        + 3 * (component_x[:, None] * ux_in + component_y[:, None] * uy_in)
        + 9 * (component_x[:, None] * ux_in + component_y[:, None] * uy_in) **2/2
        - 3 * (ux_in ** 2 + uy_in ** 2)[None, :] /2
    )


    F[:, 0, :] = F_inlet_eq
    F[:, -1, :] = F[:, -2, :]

    #F = cp.asnumpy(F)

    return F

def collision_equilibrium_CUDA(F, component_x, component_y, w, t_delta, tau):

    '''component_x = cp.asarray(component_x)
    component_y = cp.asarray(component_y)
    w = cp.asarray(w)
    F = cp.asarray(F)'''

    rho_total = cp.sum(F, axis = 0)
    ux = (F * component_x.reshape((9, 1, 1))).sum(axis=0) / rho_total
    uy = (F * component_y.reshape((9, 1, 1))).sum(axis=0) / rho_total
    u_square = ux**2 + uy**2

    Feq = np.zeros_like(F)
    for j in range(9):
        cu = component_x[j] * ux + component_y[j] * uy
        Feq[j, :, :] = w[j] * rho_total * (1 + 3 * cu
                                           + 4.5 * cu**2
                                           - 1.5 * u_square)

    F = F * (1 - t_delta / tau) + Feq * (t_delta / tau)

    #F = cp.asnumpy(F)


    return F, Feq

def collision_equilibrium_CUDA2(F, component_x, component_y, w, t_delta, tau):

    rho_total = cp.sum(F, axis = 0)
    ux = (F * component_x.reshape((9, 1, 1))).sum(axis=0) / rho_total
    uy = (F * component_y.reshape((9, 1, 1))).sum(axis=0) / rho_total
    u_square = ux**2 + uy**2

    cu = component_x[:, None, None] * ux[None, :, :] + component_y[:, None, None] * uy[None, :, :]

    Feq = w[:, None, None] * rho_total * (
        1
        + 3 * cu
        + 4.5 * cu**2
        - 1.5 * u_square
    )

    F = F * (1 - t_delta / tau) + Feq * (t_delta / tau)
    return F, Feq

def boundary_object_CUDA(F, objects):
    # this one is really slow even when it's all done in gpu
    inside_boundary = F[:, objects]
    inside_boundary = inside_boundary[[0,5,6,7,8,1,2,3,4], :]
    F[:, objects] = inside_boundary
    return F


