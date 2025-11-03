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


def compute_multispecies_moments(F, component_x, component_y, eps=1e-12):
    """
    Compute density and velocity of each species together with the barycentric velocity.

    Parameters
    ----------
    F : cupy.ndarray
        Distribution functions with shape (num_species, q, nx, ny).
    component_x, component_y : cupy.ndarray
        Lattice velocity components of length q.
    eps : float
        Small regularisation to avoid division by zero.

    Returns
    -------
    rho_sigma : cupy.ndarray
        Species densities, shape (num_species, nx, ny).
    ux_sigma, uy_sigma : cupy.ndarray
        Species velocities, same shape as rho_sigma.
    rho_total : cupy.ndarray
        Mixture density, shape (nx, ny).
    ux_mix, uy_mix : cupy.ndarray
        Barycentric velocity of the mixture, shape (nx, ny).
    """
    rho_sigma = cp.sum(F, axis=1)

    momentum_x = cp.sum(F * component_x[None, :, None, None], axis=1)
    momentum_y = cp.sum(F * component_y[None, :, None, None], axis=1)

    ux_sigma = cp.divide(momentum_x, rho_sigma, out=cp.zeros_like(momentum_x), where=rho_sigma > eps)
    uy_sigma = cp.divide(momentum_y, rho_sigma, out=cp.zeros_like(momentum_y), where=rho_sigma > eps)

    rho_total = cp.sum(rho_sigma, axis=0)

    ux_mix_num = cp.sum(rho_sigma * ux_sigma, axis=0)
    uy_mix_num = cp.sum(rho_sigma * uy_sigma, axis=0)

    ux_mix = cp.divide(ux_mix_num, rho_total, out=cp.zeros_like(ux_mix_num), where=rho_total > eps)
    uy_mix = cp.divide(uy_mix_num, rho_total, out=cp.zeros_like(uy_mix_num), where=rho_total > eps)

    return rho_sigma, ux_sigma, uy_sigma, rho_total, ux_mix, uy_mix


def equilibrium_multispecies(rho_sigma, phi, component_x, component_y, w, ux_mix, uy_mix):
    """
    Maxwellian equilibrium with species-specific compressibility factor phi.

    Parameters
    ----------
    rho_sigma : cupy.ndarray
        Species densities, shape (num_species, nx, ny).
    phi : cupy.ndarray
        Species compressibility factors, length num_species.
    component_x, component_y : cupy.ndarray
        Lattice velocity components of length q.
    w : cupy.ndarray
        Lattice weights of length q.
    ux_mix, uy_mix : cupy.ndarray
        Barycentric velocity shared by the mixture, shape (nx, ny).

    Returns
    -------
    cp.ndarray
        Equilibrium distributions with shape (num_species, q, nx, ny).
    """
    num_species, nx, ny = rho_sigma.shape
    q = component_x.shape[0]

    phi_b = phi[:, None, None]
    u_square = ux_mix**2 + uy_mix**2

    cu = cp.broadcast_to(
        component_x[None, :, None, None] * ux_mix[None, None, :, :]
        + component_y[None, :, None, None] * uy_mix[None, None, :, :]
        ,
        (num_species, q, nx, ny)
    )

    feq = cp.empty((num_species, q, nx, ny), dtype=rho_sigma.dtype)

    feq[:, 0, :, :] = (4.0 / 9.0) * rho_sigma * (
        (9.0 - 5.0 * phi_b) / 4.0 - 1.5 * u_square[None, :, :]
    )

    for idx in range(1, q):
        feq[:, idx, :, :] = w[idx] * rho_sigma * (
            phi_b
            + 3.0 * cu[:, idx, :, :]
            + 4.5 * cu[:, idx, :, :] ** 2
            - 1.5 * u_square[None, :, :]
        )

    return feq


def collision_multispecies_BGK(F, component_x, component_y, w, t_delta, tau, phi):
    """
    BGK collision for multiple species sharing a barycentric velocity.

    Parameters
    ----------
    F : cupy.ndarray
        Distribution functions with shape (num_species, q, nx, ny).
    component_x, component_y, w : cupy.ndarray
        Lattice descriptors of length q.
    t_delta : float
        Time step.
    tau : float or cupy.ndarray
        Relaxation time(s); scalar or array of length num_species.
    phi : cupy.ndarray
        Compressibility factor per species (same length as tau).

    Returns
    -------
    F_new : cupy.ndarray
        Post-collision distributions.
    feq : cupy.ndarray
        Equilibrium distributions used in the collision.
    rho_sigma : cupy.ndarray
        Species densities for convenience.
    rho_total : cupy.ndarray
        Mixture density.
    ux_mix, uy_mix : cupy.ndarray
        Barycentric velocity field.
    """
    rho_sigma, _, _, rho_total, ux_mix, uy_mix = compute_multispecies_moments(F, component_x, component_y)
    feq = equilibrium_multispecies(rho_sigma, phi, component_x, component_y, w, ux_mix, uy_mix)

    tau_array = cp.asarray(tau)
    if tau_array.ndim == 0:
        tau_array = cp.full((F.shape[0],), tau_array, dtype=F.dtype)

    omega = t_delta / tau_array[:, None, None, None]
    F_new = F * (1.0 - omega) + feq * omega

    return F_new, feq, rho_sigma, rho_total, ux_mix, uy_mix


