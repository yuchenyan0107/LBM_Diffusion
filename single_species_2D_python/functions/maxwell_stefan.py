import numpy as np


D2Q9_CX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=np.int64)
D2Q9_CY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=np.int64)
BB_OPPOSITE = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int64)
THETA = 0.5


def hydrodynamic_moments(distribution: np.ndarray) -> tuple[float, float, float]:
    """
    Compute density and velocity components for a D2Q9 population.
    """
    rho = np.sum(distribution)
    if rho <= 0.0:
        return 0.0, 0.0, 0.0

    ux = np.dot(D2Q9_CX, distribution) / rho
    uy = np.dot(D2Q9_CY, distribution) / rho
    return rho, ux, uy


def equilibrium_distribution(rho: float, phi: float, ux: float, uy: float) -> np.ndarray:
    """
    Maxwellian equilibrium for multispecies LBM (D2Q9).
    """
    u_sq = ux * ux + uy * uy

    feq = np.empty(9, dtype=np.float64)
    feq[0] = (4.0 / 9.0) * rho * ((9.0 - 5.0 * phi) / 4.0 - 1.5 * u_sq)

    cu = ux
    feq[1] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = uy
    feq[2] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -ux
    feq[3] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -uy
    feq[4] = (1.0 / 9.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)

    cu = ux + uy
    feq[5] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -ux + uy
    feq[6] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = -ux - uy
    feq[7] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)
    cu = ux - uy
    feq[8] = (1.0 / 36.0) * rho * (phi + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq)

    return feq


def mrt_matrix(lamd: float, lamn: float, lamb: float) -> np.ndarray:
    """
    Collision matrix A (velocity space) used by the θ-splitting scheme.
    Direct port of MIXLBM.f90::MRT.
    """
    lam0 = 0.0
    lam34 = 1.0

    A = np.zeros((9, 9), dtype=np.float64)

    A[0, 0] = lam0
    A[0, 1] = lam0 - lamb
    A[0, 2] = lam0 - lamb
    A[0, 3] = lam0 - lamb
    A[0, 4] = lam0 - lamb
    A[0, 5] = lam0 - 2.0 * lamb + lam34
    A[0, 6] = lam0 - 2.0 * lamb + lam34
    A[0, 7] = lam0 - 2.0 * lamb + lam34
    A[0, 8] = lam0 - 2.0 * lamb + lam34

    A[1, 1] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[1, 2] = lamb / 4.0 - lamn / 4.0
    A[1, 3] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[1, 4] = lamb / 4.0 - lamn / 4.0
    A[1, 5] = lamd / 2.0 + lamb / 2.0 - lam34
    A[1, 6] = -lamd / 2.0 + lamb / 2.0
    A[1, 7] = -lamd / 2.0 + lamb / 2.0
    A[1, 8] = lamd / 2.0 + lamb / 2.0 - lam34

    A[2, 1] = lamb / 4.0 - lamn / 4.0
    A[2, 2] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[2, 3] = lamb / 4.0 - lamn / 4.0
    A[2, 4] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[2, 5] = lamd / 2.0 + lamb / 2.0 - lam34
    A[2, 6] = lamd / 2.0 + lamb / 2.0 - lam34
    A[2, 7] = -lamd / 2.0 + lamb / 2.0
    A[2, 8] = -lamd / 2.0 + lamb / 2.0

    A[3, 1] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[3, 2] = lamb / 4.0 - lamn / 4.0
    A[3, 3] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[3, 4] = lamb / 4.0 - lamn / 4.0
    A[3, 5] = -lamd / 2.0 + lamb / 2.0
    A[3, 6] = lamd / 2.0 + lamb / 2.0 - lam34
    A[3, 7] = lamd / 2.0 + lamb / 2.0 - lam34
    A[3, 8] = -lamd / 2.0 + lamb / 2.0

    A[4, 1] = lamb / 4.0 - lamn / 4.0
    A[4, 2] = -lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[4, 3] = lamb / 4.0 - lamn / 4.0
    A[4, 4] = lamd / 2.0 + lamb / 4.0 + lamn / 4.0
    A[4, 5] = -lamd / 2.0 + lamb / 2.0
    A[4, 6] = -lamd / 2.0 + lamb / 2.0
    A[4, 7] = lamd / 2.0 + lamb / 2.0 - lam34
    A[4, 8] = lamd / 2.0 + lamb / 2.0 - lam34

    A[5, 5] = lamn / 4.0 + 3.0 / 4.0 * lam34
    A[5, 6] = -lamn / 4.0 + lam34 / 4.0
    A[5, 7] = lamn / 4.0 - lam34 / 4.0
    A[5, 8] = -lamn / 4.0 + lam34 / 4.0

    A[6, 5] = -lamn / 4.0 + lam34 / 4.0
    A[6, 6] = lamn / 4.0 + 3.0 / 4.0 * lam34
    A[6, 7] = -lamn / 4.0 + lam34 / 4.0
    A[6, 8] = lamn / 4.0 - lam34 / 4.0

    A[7, 5] = lamn / 4.0 - lam34 / 4.0
    A[7, 6] = -lamn / 4.0 + lam34 / 4.0
    A[7, 7] = lamn / 4.0 + 3.0 / 4.0 * lam34
    A[7, 8] = -lamn / 4.0 + lam34 / 4.0

    A[8, 5] = -lamn / 4.0 + lam34 / 4.0
    A[8, 6] = lamn / 4.0 - lam34 / 4.0
    A[8, 7] = -lamn / 4.0 + lam34 / 4.0
    A[8, 8] = lamn / 4.0 + 3.0 / 4.0 * lam34

    return A


def theta_transformation(lamd: float, lamn: float, lamb: float, theta: float = THETA) -> np.ndarray:
    """
    (I + theta A)^-1 matrix used in the θ-splitting solver.
    """
    lam0 = 0.0
    lam34 = 1.0

    G = np.zeros((9, 9), dtype=np.float64)

    lamd_theta = lamd * theta
    lamn_theta = lamn * theta
    lamb_theta = lamb * theta
    lam34_theta = lam34 * theta

    def denom(*vals: float) -> float:
        prod = 1.0
        for value in vals:
            prod *= (value + 1.0)
        return prod

    G[0, 0] = 1.0 / (1.0 + lam0 * theta)
    coeff_01 = -theta * (lam0 - lamb) / ((theta * lamb) + 1.0) / (1.0 + lam0 * theta)
    for idx in (1, 2, 3, 4):
        G[0, idx] = coeff_01
    coeff_05 = -(
        -theta * lamb * lam34
        + 2.0 * theta * lam0 * lam34
        - theta * lam0 * lamb
        + lam0
        - 2.0 * lamb
        + lam34
    ) * theta / ((theta * lamb) + 1.0) / ((theta * lam34) + 1.0) / (1.0 + lam0 * theta)
    for idx in (5, 6, 7, 8):
        G[0, idx] = coeff_05

    common1 = denom(lamd_theta, lamb_theta, lamn_theta) * 4.0
    G[1, 1] = (
        theta * theta * lamd * lamn
        + theta * theta * lamb * lamd
        + 2.0 * theta * theta * lamn * lamb
        + 2.0 * theta * lamd
        + 3.0 * theta * lamn
        + 3.0 * theta * lamb
        + 4.0
    ) / common1
    G[1, 2] = -(lamb - lamn) * theta / ((theta * lamn) + 1.0) / ((theta * lamb) + 1.0) / 4.0
    G[1, 3] = (
        theta * lamd * lamn
        + theta * lamb * lamd
        - 2.0 * theta * lamn * lamb
        - lamn
        - lamb
        + 2.0 * lamd
    ) * theta / ((theta * lamn) + 1.0) / ((theta * lamb) + 1.0) / (theta * lamd + 1.0) / 4.0
    G[1, 4] = G[1, 2]
    G[1, 5] = (
        theta * lamd * lam34
        - 2.0 * theta * lamb * lamd
        + theta * lamb * lam34
        - lamd
        + 2.0 * lam34
        - lamb
    ) * theta / ((theta * lamb) + 1.0) / ((theta * lam34) + 1.0) / (theta * lamd + 1.0) / 2.0
    G[1, 6] = (lamd - lamb) * theta / ((theta * lamb) + 1.0) / (theta * lamd + 1.0) / 2.0
    G[1, 7] = G[1, 6]
    G[1, 8] = G[1, 5]

    G[2, 1] = G[1, 2]
    G[2, 2] = G[1, 1]
    G[2, 3] = G[1, 2]
    G[2, 4] = G[1, 3]
    G[2, 5] = G[1, 5]
    G[2, 6] = G[1, 5]
    G[2, 7] = G[1, 6]
    G[2, 8] = G[1, 6]

    G[3, 1] = G[1, 3]
    G[3, 2] = G[1, 2]
    G[3, 3] = G[1, 1]
    G[3, 4] = G[1, 2]
    G[3, 5] = G[1, 6]
    G[3, 6] = G[1, 5]
    G[3, 7] = G[1, 5]
    G[3, 8] = G[1, 6]

    G[4, 1] = G[1, 2]
    G[4, 2] = G[1, 3]
    G[4, 3] = G[1, 2]
    G[4, 4] = G[1, 1]
    G[4, 5] = G[1, 6]
    G[4, 6] = G[1, 6]
    G[4, 7] = G[1, 5]
    G[4, 8] = G[1, 5]

    base = ((theta * lamn) + 1.0) * ((theta * lam34) + 1.0) * 4.0
    G[5, 5] = (theta * lam34 + 3.0 * theta * lamn + 4.0) / base
    G[5, 6] = -(lam34 - lamn) * theta / base
    G[5, 7] = (lam34 - lamn) * theta / base
    G[5, 8] = G[5, 6]

    G[6, 5] = -G[5, 6]
    G[6, 6] = G[5, 5]
    G[6, 7] = -G[5, 6]
    G[6, 8] = G[5, 7]

    G[7, 5] = (lam34 - lamn) * theta / base
    G[7, 6] = -G[5, 6]
    G[7, 7] = G[5, 5]
    G[7, 8] = -G[5, 6]

    G[8, 5] = -G[5, 6]
    G[8, 6] = G[5, 7]
    G[8, 7] = -G[5, 6]
    G[8, 8] = G[5, 5]

    return G


def B_binary_resistivity(m1: float, m2: float, nB: int) -> float:
    table = np.array(
        [
            0.500000000000000,
            0.601274801898418,
            0.723062774795962,
            0.869518853351124,
            1.045639552591273,
            1.257433429682934,
            1.512126072666108,
            1.818406609575491,
            2.186724147886553,
            2.629644257653945,
            3.162277660168372,
            3.802795747331057,
            4.573050519273251,
            5.499320090094956,
            6.613205195495662,
            7.952707287670479,
            9.563524997900334,
            11.500613197126169,
            13.830057843624722,
            16.631330580338215,
        ],
        dtype=np.float64,
    )
    factor = table[nB - 1]
    return factor * 10.0 * (1.0 / m1 + 1.0 / m2) ** (-0.5)


def D_fick_diffusivity(m: float) -> float:
    table = np.array(
        [
            0.010000000000000,
            0.013645831365889,
            0.018620871366629,
            0.025409727055493,
            0.034673685045253,
            0.047315125896148,
            0.064565422903466,
            0.088104887300801,
            0.120226443461741,
            0.164058977319954,
            0.223872113856834,
            0.305492111321552,
            0.416869383470336,
            0.568852930843842,
            0.776247116628692,
            1.059253725177290,
            1.445439770745929,
            1.972422736114855,
            2.691534803926918,
            3.672823004980850,
        ],
        dtype=np.float64,
    )
    return table[4] * 0.2 / m


def NU_kinematic_viscosity(nN: int) -> float:
    table = np.array(
        [
            1.0000,
            1.1708,
            1.3707,
            1.6048,
            1.8789,
            2.1998,
            2.5754,
            3.0153,
            3.5302,
            4.1331,
            4.8390,
            5.6654,
            6.6329,
            7.7657,
            9.0919,
            10.6446,
            12.4625,
            14.5908,
            17.0826,
            20.0000,
        ],
        dtype=np.float64,
    )
    return table[nN - 1]


def XI_bulk_viscosity() -> float:
    return 0.4


def solve_linear_system(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    return np.linalg.solve(matrix, rhs)


class MaxwellStefanState:
    """
    Container for the lattice populations and macroscopic fields.
    """

    def __init__(self, num_species: int, nx: int, ny: int):
        self.fd = np.zeros((num_species, 9, nx, ny), dtype=np.float64)
        self.fd_new = np.zeros_like(self.fd)
        self.rsigma = np.zeros((num_species, nx, ny), dtype=np.float64)
        self.psigma = np.zeros_like(self.rsigma)
        self.uxsigma = np.zeros_like(self.rsigma)
        self.uysigma = np.zeros_like(self.rsigma)


def initialise_from_partial_pressures(
    state: MaxwellStefanState,
    phi: np.ndarray,
    psigma_init: np.ndarray,
) -> None:
    num_species, nx, ny = psigma_init.shape
    for s in range(num_species):
        rho_sigma = 3.0 * psigma_init[s] / phi[s]
        feq = equilibrium_distribution(1.0, phi[s], 0.0, 0.0)
        for i in range(nx):
            for j in range(ny):
                feq_local = equilibrium_distribution(rho_sigma[i, j], phi[s], 0.0, 0.0)
                state.fd[s, :, i, j] = feq_local
                state.fd_new[s, :, i, j] = feq_local
                state.rsigma[s, i, j], state.uxsigma[s, i, j], state.uysigma[s, i, j] = hydrodynamic_moments(
                    feq_local
                )
                state.psigma[s, i, j] = phi[s] * state.rsigma[s, i, j] / 3.0


def maxwell_stefan_step(
    state: MaxwellStefanState,
    phi: np.ndarray,
    molecular_weights: np.ndarray,
    nN: int = 20,
    nB: int = 15,
    theta: float = THETA,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform one θ-splitting Maxwell–Stefan collision/streaming step.
    This is a direct translation of MIXLBM.f90 with periodic boundaries.
    """
    num_species, _, nx, ny = state.fd.shape

    fd = state.fd
    fd_new = np.zeros_like(fd)

    rsigma_local = np.zeros(num_species, dtype=np.float64)
    psigma_local = np.zeros(num_species, dtype=np.float64)
    uxsigma_local = np.zeros(num_species, dtype=np.float64)
    uysigma_local = np.zeros(num_species, dtype=np.float64)
    uxstar = np.zeros(num_species, dtype=np.float64)
    uystar = np.zeros(num_species, dtype=np.float64)
    lambda_sigma = np.zeros(num_species, dtype=np.float64)
    lamd_sigma = np.zeros(num_species, dtype=np.float64)
    lamn_sigma = np.zeros(num_species, dtype=np.float64)
    lamb_sigma = np.zeros(num_species, dtype=np.float64)
    CHI = np.zeros((num_species, num_species), dtype=np.float64)

    nu_value = NU_kinematic_viscosity(nN)
    xi_value = XI_bulk_viscosity()

    for i in range(nx):
        for j in range(ny):
            for k in range(9):
                iI = (i + D2Q9_CX[k]) % nx
                jI = (j + D2Q9_CY[k]) % ny

                for s in range(num_species):
                    f_neighbor = fd[s, :, iI, jI]
                    rho, ux, uy = hydrodynamic_moments(f_neighbor)
                    rsigma_local[s] = rho
                    psigma_local[s] = phi[s] * rho / 3.0
                    uxsigma_local[s] = ux
                    uysigma_local[s] = uy

                r_mix = np.sum(rsigma_local)
                p_mix = np.sum(psigma_local)
                if r_mix <= 0.0:
                    r_mix = 1.0
                if p_mix <= 0.0:
                    p_mix = 1.0

                ux_mix = np.dot(rsigma_local, uxsigma_local) / r_mix
                uy_mix = np.dot(rsigma_local, uysigma_local) / r_mix

                inv_MM_sum = 0.0
                for s in range(num_species):
                    inv_MM_sum += (rsigma_local[s] / r_mix) / molecular_weights[s]
                mixture_mass = 1.0 / inv_MM_sum if inv_MM_sum != 0.0 else 0.0

                for s in range(num_species):
                    for vs in range(num_species):
                        resistor = B_binary_resistivity(
                            molecular_weights[s],
                            molecular_weights[vs],
                            nB,
                        )
                        ref = B_binary_resistivity(
                            molecular_weights[s],
                            molecular_weights[s],
                            nB,
                        )
                        CHI[s, vs] = (
                            (mixture_mass * mixture_mass)
                            / (molecular_weights[s] * molecular_weights[vs])
                            * resistor
                            / ref
                        )

                for s in range(num_species):
                    uxstar[s] = uxsigma_local[s]
                    uystar[s] = uysigma_local[s]
                    for vs in range(num_species):
                        factor = CHI[s, vs] * (rsigma_local[vs] / r_mix)
                        uxstar[s] += factor * (uxsigma_local[vs] - uxsigma_local[s])
                        uystar[s] += factor * (uysigma_local[vs] - uysigma_local[s])

                for s in range(num_species):
                    lambda_sigma[s] = p_mix * B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[s],
                        nB,
                    ) / r_mix
                    lamd_sigma[s] = lambda_sigma[s]
                    lamn_sigma[s] = 1.0 / (3.0 * nu_value)
                    lamb_sigma[s] = (2.0 - phi[s]) / (3.0 * xi_value)

                for s in range(num_species):
                    f_neighbor = fd[s, :, iI, jI]
                    feq = equilibrium_distribution(
                        rsigma_local[s],
                        phi[s],
                        uxstar[s],
                        uystar[s],
                    )
                    A = mrt_matrix(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s])
                    delta = feq - f_neighbor
                    f_star = f_neighbor - theta * (A @ delta)
                    delta_star = feq - f_star
                    f_coll = f_star + A @ delta_star
                    fd_new[s, BB_OPPOSITE[k], i, j] = f_coll[BB_OPPOSITE[k]]

    rsigma_post = np.zeros((num_species, nx, ny), dtype=np.float64)
    uxsigma_post = np.zeros_like(rsigma_post)
    uysigma_post = np.zeros_like(rsigma_post)
    psigma_post = np.zeros_like(rsigma_post)

    for i in range(nx):
        for j in range(ny):
            for s in range(num_species):
                rho, ux, uy = hydrodynamic_moments(fd_new[s, :, i, j])
                rsigma_post[s, i, j] = rho
                uxsigma_post[s, i, j] = ux
                uysigma_post[s, i, j] = uy
                psigma_post[s, i, j] = phi[s] * rho / 3.0

            r_mix = np.sum(rsigma_post[:, i, j])
            p_mix = np.sum(psigma_post[:, i, j])
            if r_mix <= 0.0:
                continue
            if p_mix <= 0.0:
                p_mix = 1.0

            inv_MM_sum = 0.0
            for s in range(num_species):
                inv_MM_sum += (rsigma_post[s, i, j] / r_mix) / molecular_weights[s]
            mixture_mass = 1.0 / inv_MM_sum if inv_MM_sum != 0.0 else 0.0

            for s in range(num_species):
                lambda_sigma[s] = p_mix * B_binary_resistivity(
                    molecular_weights[s],
                    molecular_weights[s],
                    nB,
                ) / r_mix
                lamd_sigma[s] = lambda_sigma[s]
                lamn_sigma[s] = 1.0 / (3.0 * nu_value)
                lamb_sigma[s] = (2.0 - phi[s]) / (3.0 * xi_value)

            for s in range(num_species):
                for vs in range(num_species):
                    resistor = B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[vs],
                        nB,
                    )
                    ref = B_binary_resistivity(
                        molecular_weights[s],
                        molecular_weights[s],
                        nB,
                    )
                    CHI[s, vs] = (
                        (mixture_mass * mixture_mass)
                        / (molecular_weights[s] * molecular_weights[vs])
                        * resistor
                        / ref
                    )

            CHIsigma = np.zeros(num_species, dtype=np.float64)
            for s in range(num_species):
                temp = 0.0
                for vs in range(num_species):
                    temp += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix)
                CHIsigma[s] = temp

            A_matrix = np.zeros((num_species, num_species), dtype=np.float64)
            for s in range(num_species):
                for vs in range(num_species):
                    if s == vs:
                        A_matrix[s, vs] = 1.0 + theta * lambda_sigma[s] * CHIsigma[s]
                    else:
                        A_matrix[s, vs] = 0.0
                    A_matrix[s, vs] -= (
                        theta
                        * lambda_sigma[s]
                        * (rsigma_post[s, i, j] / r_mix)
                        * CHI[s, vs]
                    )

            gjx = rsigma_post[:, i, j] * uxsigma_post[:, i, j]
            gjy = rsigma_post[:, i, j] * uysigma_post[:, i, j]

            jx = solve_linear_system(A_matrix, gjx)
            jy = solve_linear_system(A_matrix, gjy)

            for s in range(num_species):
                uxsigma_post[s, i, j] = jx[s] / rsigma_post[s, i, j] if rsigma_post[s, i, j] > 0 else 0.0
                uysigma_post[s, i, j] = jy[s] / rsigma_post[s, i, j] if rsigma_post[s, i, j] > 0 else 0.0

            for s in range(num_species):
                uxstar[s] = uxsigma_post[s, i, j]
                uystar[s] = uysigma_post[s, i, j]
                for vs in range(num_species):
                    uxstar[s] += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix) * (
                        uxsigma_post[vs, i, j] - uxsigma_post[s, i, j]
                    )
                    uystar[s] += CHI[s, vs] * (rsigma_post[vs, i, j] / r_mix) * (
                        uysigma_post[vs, i, j] - uysigma_post[s, i, j]
                    )

            for s in range(num_species):
                feq = equilibrium_distribution(rsigma_post[s, i, j], phi[s], uxstar[s], uystar[s])
                A = mrt_matrix(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s])
                G = theta_transformation(lamd_sigma[s], lamn_sigma[s], lamb_sigma[s], theta=theta)
                tempV = theta * (A @ feq)
                tempV1 = fd_new[s, :, i, j] + tempV
                fd_new[s, :, i, j] = G @ tempV1

    state.fd = fd_new
    state.fd_new = fd_new.copy()
    state.rsigma = rsigma_post
    state.psigma = psigma_post
    state.uxsigma = uxsigma_post
    state.uysigma = uysigma_post

    rho_mix = np.sum(rsigma_post, axis=0)
    ux_mix = np.sum(rsigma_post * uxsigma_post, axis=0) / np.where(rho_mix > 0.0, rho_mix, 1.0)
    uy_mix = np.sum(rsigma_post * uysigma_post, axis=0) / np.where(rho_mix > 0.0, rho_mix, 1.0)

    return rsigma_post, ux_mix, uy_mix
