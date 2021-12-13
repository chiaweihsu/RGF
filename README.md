# RGF
RGF is a finite-difference frequency-domain (FDFD) solver that uses the <ins>r</ins>ecursive <ins>G</ins>reen's <ins>f</ins>unction method to solve the scattering matrix (or a portion of it) of a user-specified structure. It solves Maxwell's equations for 2D TM modes (`cal_smatrix_RGF.m`) or the 3D scalar wave equation (`cal_smatrix_RGF_3D.m`). The user-specified structure can have an arbitrary permittivity profile (which can be complex-valued to describe linear gain or absorption) and is surrounded by semi-infinite homogeneous spaces on the left and right sides (with user-specified constant permittivity) where the input and output channels are defined. In the transverse direction(s), the user can specify periodic, Bloch periodic, or Dirichlet boundary conditions .

RGF solves the discretized wave equation with no additional approximation. The scattering matrix is normalized by the longitudinal flux, so that the full scattering matrix is unitary if the refractive index profile is real-valued. RGF supports complex-valued frequency, i.e. it can return the scattering matrix analytically continued from a real frequency to an arbitrary complex frequency.

Detailed usage and explanation are given in the first contiguous comment lines, also shown with command <code>help cal_smatrix_RGF</code> or <code>help cal_smatrix_RGF_3D</code> in MATLAB. For basic usage, see `example_script.m` and `example_script_3D.m`.

[Michael Wimmer's PhD thesis](https://epub.uni-regensburg.de/12142/) is a good reference on the recursive Green's function method.

You can cite this software as:
C. W. Hsu, RGF, https://github.com/chiaweihsu/RGF
