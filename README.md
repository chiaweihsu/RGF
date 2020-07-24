# RGF
RGF is a finite-difference frequency-domain (FDFD) solver that uses the Recursive Green's function method to solve the scattering matrix (or a portion of it) of a user-specified structure. It solves Maxwell's equations for 2D TM modes discretized with finite difference. The user-specified structure can have an arbitrary refractive index profile (which can be complex-valued to describe gain or absorption) and is surrounded by infinite homogeneous spaces on the left and right sides, with either periodic or Dirichlet boundary condition above and below.

RGF solves the discretized wave equation with no additional approximation. The scattering matrix is normalized by the longitudinal flux, so that the full scattering matrix is unitary if the refractive index profile is real-valued. RGF supports complex-valued frequency, i.e. it can return the scattering matrix analytically continued from a real frequency to an arbitrary complex frequency.

Detailed usage and explanation are given in the first contiguous comment lines in cal_smatrix_RGF.m, also shown with command "help cal_smatrix_RGF" in MATLAB. For common usage, see example_script.m.
