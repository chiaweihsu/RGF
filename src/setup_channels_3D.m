function channels = setup_channels_3D(nx, ny, BC, k0dx, epsilon_L, epsilon_R)
%SETUP_CHANNELS_3D    Set up channels of the homogeneous spaces.
%   channels = setup_channels_3D(nx, ny, BC, k0dx, epsilon_L, epsilon_R) returns a
%   structure containing properties of the propagating and evanescent channels
%   in the homogeneous spaces on the left and right.
%
%   === Input parameters ===
%   nx (positive integer scalar):
%       Number of grid points in the x direction.
%   ny (positive integer scalar):
%       Number of grid points in the y direction.
%   BC (character vector or string (case insensitive) or real vector):
%       Boundary condition in the x and y directions.
%       For character inputs, the options are:
%         'periodic': f(n+nx,m) = f(n,m+ny) = f(n,m)
%         'Dirichlet' or 'PEC': f(0,m) = f(nx+1,m) = f(n,0) = f(n,ny+1) = 0
%       When BC is a real vector, the Bloch periodic boundary condition is used
%       with f(n+nx,m) = f(n,m)*exp(1i*kxa), f(n,m+ny) = f(n,m)*exp(1i*kya).
%       Here, BC = [kxa, kya]
%   k0dx (numeric scalar, real or complex):
%       Normalized frequency k0*dx = (2*pi/lambda)*dx.
%   epsilon_L (numeric scalar, real or complex):
%       Relative permittivity of the homogeneous space on the left.
%   epsilon_R (numeric scalar, real or complex, optional):
%       Relative permittivity of the homogeneous space on the right. Only the
%       left side will be considered if epsilon_R is not given or if it is empty
%       or NaN.
%
%   === Returns ===
%   channels (structure):
%       channels.kxdx (nx-by-ny real matrix):
%           Normalized transverse wave number kx*dx for all nx*ny channels. They
%           are real-valued and are ordered from small to large.
%       channels.kydx (nx-by-ny real matrix):
%           Normalized transverse wave number ky*dx for all nx*ny channels. They
%           are real-valued and are ordered from small to large.
%       channels.fun_phi (function_handle):
%           A function that, given one element of kxdx and kydx as the input,
%           returns its normalized transverse field profile as an nx-by-ny matrix.
%       channels.L and channels.R (structue):
%           Structures containing properties specific to the left (L) and right
%           (R) sides.
%       channels.L.N_prop (integer scalar):
%           Number of propagating channels.
%       channels.L.kzdx (nx-by-ny complex matrix):
%           Normalized longitudinal wave number kz*dx for all channels. Due to
%           the discretization, kzdx is equivalent to kzdx + 2*pi, so kzdx is
%           on the folded complex plane (ie, the surface of a cylinder).
%           There are 2*nx*ny unique solutions of kzdx, and whenever kzdx is a
%           solution, -kzdx is also a solution. Here, we choose the nx*ny solutions
%           such that
%           1) When k0dx is real, we have
%               Propagating channels: 0 < Re(kzdx) < pi, Im(kzdx) = 0.
%               Evanescent channels: Im(kzdx) >= 0, mod(Re(kzdx),pi) = 0.
%           2) When k0dx is complex, we analytically continue the above choice
%           onto the complex-frequency plane. Specifically, we pick the kzdx
%           that is continuously connected to the one with real k0dx through a
%           vertical line in the complex-(k0dx^2) plane.
%       channels.L.sqrt_mu (1-by-N_prop row vector):
%           Square root of the normalized longitudinal group velocity of the
%           propagating channels, sqrt_mu = sqrt(sin(kzdx)). The longitudinal
%           group velocity is v_g = (sin(kzdx)/k0dx)*(c/epsilon_L).
%       channels.L.ind_prop (1-by-N_prop integer row vector):
%           Indices of the N_prop propagating channels among all ny channels.

% Check input parameters
if ~(isscalar(nx) && isnumeric(nx) && isreal(nx) && (round(nx)==nx) && nx>0); error('nx must be a positive integer scalar'); end
if ~(isscalar(ny) && isnumeric(ny) && isreal(ny) && (round(ny)==ny) && ny>0); error('ny must be a positive integer scalar'); end
if ~((ischar(BC) && isvector(BC)) || (isstring(BC) && isscalar(BC)) || (isnumeric(BC) && isvector(BC) && numel(BC)==2))
    error('Boundary condition BC has to be a character vector or string, or numeric vector with 2 elements (for Bloch periodic BC)');
end
if ~(isscalar(k0dx) && isnumeric(k0dx)); error('k0dx must be a numeric scalar'); end
if ~(isscalar(epsilon_L) && isnumeric(epsilon_L)); error('epsilon_L must be a numeric scalar'); end
if nargin == 5 || isempty(epsilon_R) || isnan(epsilon_R)
    two_sided = false;
else
    if ~(isscalar(epsilon_R) && isnumeric(epsilon_R)); error('epsilon_R must be a numeric scalar'); end
    two_sided = true;
end

% handle periodic and Bloch periodic boundary conditions
if strcmpi(BC, 'Bloch')
    error('To use Bloch periodic boundary condition, set BC to [kx*ax, ky*ay] where (kx,ky) is the Bloch wave vector and (ax,ay) is the periodicity');
elseif isnumeric(BC)
    kxa = BC(1);
    kya = BC(2);
    BC = 'Bloch';
    % kxa and kya must be real for matrix phi to be unitary
    if ~isreal(kxa)
        error('kx*a = %g + 1i*%g is a complex number; has to be real for a complete orthonormal transverse basis.', real(kxa), imag(kxa));
    end
    if ~isreal(kya)
        error('ky*a = %g + 1i*%g is a complex number; has to be real for a complete orthonormal transverse basis.', real(kya), imag(kya));
    end
elseif strcmpi(BC, 'periodic')
    kxa = 0;
    kya = 0;
    BC = 'Bloch';
end

% Transverse modes (form a complete basis and are independent of epsilon_L/R)
if strcmpi(BC, 'Bloch')
    % The transverse mode index where kxdx = kxa/ny
    if mod(nx,2) == 1
        ind_zero_kx = round((nx+1)/2);
    else
        ind_zero_kx = round(nx/2);
    end
    % The transverse mode index where kydx = kya/ny
    if mod(ny,2) == 1
        ind_zero_ky = round((ny+1)/2);
    else
        ind_zero_ky = round(ny/2);
    end
    channels.kxdx = repmat((kxa/nx) + ((1:nx).'-ind_zero_kx)*(2*pi/nx), 1, ny);
    channels.kydx = repmat((kya/ny) + ((1:ny)-ind_zero_ky)*(2*pi/ny), nx, 1);
    % Normalized transverse mode profile; use implicit expansion, and let x be the fast-changing index
    channels.fun_phi = @(kxdx,kydx) exp(((1:nx).')*(1i*kxdx)).*(exp(((1:ny))*(1i*kydx))/sqrt(nx*ny));
elseif strcmpi(BC, 'Dirichlet') || strcmpi(BC, 'PEC') % PEC on both sides
    channels.kxdx = repmat(((1:nx).')*(pi/(nx+1)), 1, ny);
    channels.kydx = repmat((1:ny)*(pi/(ny+1)), nx, 1);
    % Normalized transverse mode profile; use implicit expansion, and let x be the fast-changing index
    channels.fun_phi = @(kxdx,kydx) sin(((1:nx).')*kxdx).*(sin(((1:ny))*kydx)*sqrt(4/((nx+1)*(ny+1)))); 
else
    error('Boundary condition BC = %s is not an eligible option', BC);
end

% Longitudinal properties for homogeneous space on the left (kzdx, sqrt_mu, number of propagating channels, etc; depends on epsilon_L/R)
channels.L = setup_longitudinal((k0dx^2)*epsilon_L, channels.kxdx, channels.kydx);

% Homogeneous space on the right
if two_sided
    if epsilon_R == epsilon_L
        channels.R = channels.L;
    else
        channels.R = setup_longitudinal((k0dx^2)*epsilon_R, channels.kxdx, channels.kydx);
    end
elseif nargin == 6
    % create channels.R.N_prop as this variable may be accessed even in one-sided systems
    channels.R.N_prop = NaN;
end

end


function side = setup_longitudinal(k0dx2_epsilon, kxdx, kydx)
% Returns a structure 'side'. See comments at the beginning of this file for more info.

% finite-difference dispersion for homogeneous space: k0dx2_epsilon = 4*sin^2(kxdx/2) + 4*sin^2(kydx/2) + 4*sin^2(kzdx/2)

% sin_kzdx_over_two_sq = sin^2(kzdx/2)
sin_kzdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kxdx/2).^2 - sin(kydx/2).^2;

% Normalized longitudinal wave number
% asin(sqrt(z)) has two branch points (at z=0 and z=1) and with the branch cuts going outward on the real-z axis; we will address the branch choice below
% Note kzdx is only defined up to modulo 2*pi (ie, kzdx is equivalent to kzdx + 2*pi, kzdx - 2*pi, etc) because sin(kzdx) and exp(1i*kzdx) are both invariant under 2pi shifts.
side.kzdx = 2*asin(sqrt(sin_kzdx_over_two_sq));

% Indices of the propagating channels
% When k0dx2_epsilon is real, these are indicies of the channels with real-valued kzdx
% When k0dx2_epsilon is complex, these are indicies of the channels we consider "propagating-like"; they have complex kzdx with 0 < real(kzdx) < pi. When k0dx2_epsilon is tuned to a real number continuously, this set continuously becomes that at real k0dx2_epsilon.
side.ind_prop = reshape(find((real(sin_kzdx_over_two_sq) > 0) & (real(sin_kzdx_over_two_sq) < 1)), 1, []);

% Number of propagating channels
side.N_prop = length(side.ind_prop);
if side.N_prop==0 && length(kydx)==1; side.ind_prop = zeros(1,0); end  % a rare possibility, but in this case ind_prop would be zeros(0,0) while it should be zeros(1,0)

% Here we address the sign choice of kzdx, namely its branch
% When k0dx2_epsilon is real, we choose the sign of kzdx such that:
% 1) 0 < kzdx < pi for propagating channels (where kzdx is real)
% 2) Im(kzdx) >= 0 for evanescent channels
% Using the correct sign is important when we build the retarded Green's function of the semi-infinite homogeneous space.
% The default branch choice of asin(sqrt(z)) returns the correct sign for the most part, except when z > 1. We need to flip the sign of those (which can only occur if k0dx2_epsilon > 4).
% When k0dx2_epsilon is complex-valued, it is not always possible to unambiguously choose the sign that is "physical", because kzdx will be complex-valued, and the sign we "want" for real(kzdx) and the sign we want for imag(kzdx) may be incompatible.
% What we do with complex k0dx2_epsilon is that we choose the sign for the (complex-valued) kzdx such that when k0dx2_epsilon is tuned to a real number continuously by fixing Re(k0dx2_epsilon) and varying Im(k0dx2_epsilon), the kzdx we choose continuously becomes the "correct" one at real k0dx2_epsilon without crossing any branch cut. To do so, we rotate the two branch cuts of asin(sqrt(z)) by 90 degrees to the lower part of the complex-z plane (ie, the lower part of the complex-k0dx2_epsilon plane), and we pick the branch with the correct sign when k0dx2_epsilon is real. This is implemented by flipping the sign of kzdx for the ones that require flipping.
% Note that we will get a discontinuity whenever k0dx2_epsilon crosses one of those vertical-pointing branch cuts. That is unavoidable.
% The following few lines implement the "flipping".
if ~isreal(k0dx2_epsilon) || (isreal(k0dx2_epsilon) && k0dx2_epsilon > 4)
    % Note that when imag(sin_kzdx_over_two_sq)=0, flipping is needed for sin_kzdx_over_two_sq>1 but not needed for sin_kzdx_over_two_sq<0.
    ind_flip = find((real(sin_kzdx_over_two_sq)<0 & imag(sin_kzdx_over_two_sq)<0) | (real(sin_kzdx_over_two_sq)>1 & imag(sin_kzdx_over_two_sq)<=0));
    side.kzdx(ind_flip) = -side.kzdx(ind_flip);
end

% Square root of the normalized longitudinal group velocity, sqrt(sin(kzdx)), for the propagating channels
% When k0dx2_epsilon is real, sqrt_mu is also real. When k0dx2_epsilon is complex, sqrt_mu is also complex.
side.sqrt_mu = sqrt(sin(side.kzdx(side.ind_prop)));

end
