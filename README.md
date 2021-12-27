# RGF
**RGF** is a finite-difference frequency-domain (FDFD) code in MATLAB that uses the <ins>r</ins>ecursive <ins>G</ins>reen's <ins>f</ins>unction method to solve the scattering matrix (or a portion of it) of a user-specified structure. Specifically, it solves the scalar Helmholtz equation <img src="https://render.githubusercontent.com/render/math?math=[\nabla^2%2Bk^2\varepsilon_{\rm r}({\bf r})]\psi({\bf r})=0"> in 2D (`cal_smatrix_RGF.m`) or 3D (`cal_smatrix_RGF_3D.m`). The 2D one is equivalent to Maxwell's equations for 2D transverse-magnetic waves with <img src="https://render.githubusercontent.com/render/math?math=\psi(x,y)=E_z(x,y)">, and the 3D one can be used as a scalar approximation of the vectorial Maxwell's equations.

The user-specified relative permittivity profile <img src="https://render.githubusercontent.com/render/math?math=\varepsilon_{\rm r}({\bf r})"> in the scattering region can be an arbitrary real- or complex-valued matrix in 2D or 3D array in 3D, with the imaginary part describing linear gain or absorption. The scattering region is surrounded by semi-infinite homogeneous spaces on the left and right sides (with user-specified constants <img src="https://render.githubusercontent.com/render/math?math=\varepsilon_{\rm r}^{\rm L}"> and <img src="https://render.githubusercontent.com/render/math?math=\varepsilon_{\rm r}^{\rm R}">) where the input and output channels are defined. In the transverse direction(s), the user can specify periodic, Bloch periodic, or Dirichlet boundary conditions .

**RGF** solves the discretized wave equation with no additional approximation. The outgoing boundary condition is implemented exactly through the retarded Green's function of a semi-infinite space on a square lattice. The scattering matrix **RGF** computes is normalized by the longitudinal flux and is exactly unitary when <img src="https://render.githubusercontent.com/render/math?math=\varepsilon_{\rm r}({\bf r})"> is real-valued.

**RGF** supports complex-valued frequency, i.e. it can return the scattering matrix analytically continued from a real frequency to an arbitrary complex frequency.

Detailed usage and explanation are given in the first contiguous comment lines, also shown with command <code>help cal_smatrix_RGF</code> or <code>help cal_smatrix_RGF_3D</code> in MATLAB. For basic usage, see `example_script.m` and `example_script_3D.m`.

The main computations in the recursive Green's function method are multiplicatins and inversions of <img src="https://render.githubusercontent.com/render/math?math=N\times N"> dense matrices, where <img src="https://render.githubusercontent.com/render/math?math=N=W/\Delta x"> in 2D and <img src="https://render.githubusercontent.com/render/math?math=N=(W/\Delta x)^2"> for square cross sections in 3D, *W* is the transverse width, and *Δx* is the discretization grid size. So the computing time scales as <img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(LW^3)"> in 2D and <img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(LW^6)"> in 3D, where *L* is the longitudinal length. The memory usage scales as <img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(LW%2BW^2)"> in 2D and <img src="https://render.githubusercontent.com/render/math?math=\mathcal{O}(LW^2%2BW^4)"> in 3D. [Michael Wimmer's PhD thesis](https://epub.uni-regensburg.de/12142/) is a good reference on the recursive Green's function method.

You can cite this software as:
C. W. Hsu, RGF, https://github.com/chiaweihsu/RGF
