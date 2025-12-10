# LBM_Diffusion
A Lattice-Boltzmann Method tool to simulate diffusion-advection of multi-species systems.

The algorithm is largely inspired by Musubi, an open-source LBM project.
https://geb.inf.tu-dresden.de/doxy/musubi/page/features/multispecies.html

Zudrop, J., Masilamani, K., Roller, S., & Asinari, P. (2017). A robust lattice Boltzmann method for parallel simulations of multicomponent flows in complex geometries. Computers & Fluids, 153, 20–33. https://doi.org/10.1016/j.compfluid.2017.04.021

The implementation of the Dirichlet and Robin boundary conditions are based on:
Zhang, T., Shi, B., Guo, Z., Chai, Z., & Lu, J. (2012). General bounce-back scheme for concentration boundary condition in the lattice-Boltzmann method. Physical Review E, 85(1), 016701. https://doi.org/10.1103/physreve.85.016701

Here we give some examples in jupyter notebooks to demonstrate the current implementation:

1. (Almost) pure diffusion in LBM: 
we compared the result from the pure diffusion model ($\nabla^2c=0$) which people used to study selective area growth (SAG) in MOCVD, to the diffusion in LBM.
The shape and the thickness of these two models matches in general, but with some differences.
With this demo, we showed that the Dirichlet and Robin boundary conditions are correctly modeled

2. The laminar shear wave: we used the decay in flow speed to calculate the kinematic viscosity, and confirmed that it matches the theoretical value.
$\nu = c_{s}^{2}(\frac{1}{rel}-\frac{1}{2})$, where the effective relaxation $rel = \frac{\lambda}{1+\theta\lambda}$.

3. Flow around cylinder: we tried to demonstrate the turbulence phenomena with the flow around cylinder. 
The result (average flow field and Reynolds stresses) was compared with the literature. 
Nguyen, Q. D., & Lei, C. (2022). A PIV study of blockage ratio effects on flow over a confined circular cylinder at low Reynolds numbers. Experiments in Fluids, 64(1). https://doi.org/10.1007/s00348-022-03548-w
(Note: the results didn't exactly match, we are still looking into it)

4. A simple case of the diffusion of a narrow stripe. We verified the diffusion follows the square root of time.

5. Tri-species diffusion demonstrating the Maxwell-Stefan effect. 
A chamber is initialized with uniform 50% nitrogen, 50% carbon dioxide in the left half, and 50% hydrogen in the right half.
The result shows that despite the nitrogen concentration is driven by the diffusion of other species, because of differences in the Maxwell-Stefan coefficients.
(We will try to find method to validate this result quantitatively)
