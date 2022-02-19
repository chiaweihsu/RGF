function channels = setup_channels(ny, BC, k0dx, epsilon_L, epsilon_R)
%SETUP_CHANNELS    Set up channels of the homogeneous spaces.
%   channels = setup_channels(ny, BC, k0dx, epsilon_L, epsilon_R) returns a
%   structure containing properties of the propagating and evanescent channels
%   in the homogeneous spaces on the left and right.
%
%   === Input parameters ===
%   ny (positive integer scalar):
%       Number of grid points in the transverse (y) direction.
%   BC (character vector or string (case insensitive) or scalar number):
%       Boundary condition in the y direction.
%       For character inputs, the options are:
%         'periodic': f(m+ny) = f(m)
%         'Dirichlet' or 'PEC': f(0) = f(ny+1) = 0
%         'Neumann' or 'PMC': f(0) = f(1), f(ny+1) = f(ny)
%         'DirichletNeumann' or 'PECPMC': f(0) = 0, f(ny+1) = f(ny)
%         'NeumannDirichlet' or 'PMCPEC': f(0) = f(1), f(ny+1) = 0
%       When BC is a scalar number, the Bloch periodic boundary
%       condition is used with f(m+ny) = f(m)*exp(1i*ky*a) where ky
%       is the Bloch wave number. Here, BC = ky*a = ky*dx*ny.
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
%       channels.kydx (1-by-ny real row vector):
%           Normalized transverse wave number ky*dx for all ny channels. They
%           are real-valued and are ordered from small to large.
%       channels.fun_phi (function_handle):
%           A function that, given one element of kydx as the input, returns its
%           normalized transverse field profile as an ny-by-1 column vector.
%           When the input kydx is a row vector, it returns a matrix where each
%           column is the respective transverse profile. The transverse modes
%           form a complete and orthonormal set, so the ny-by-ny matrix
%           channels.fun_phi(channels.kydx) is unitary.
%       channels.L and channels.R (structue):
%           Structures containing properties specific to the left (L) and right
%           (R) sides.
%       channels.L.N_prop (integer scalar):
%           Number of propagating channels.
%       channels.L.kxdx (1-by-ny complex row vector):
%           Normalized longitudinal wave number kx*dx for all channels. Due to
%           the discretization, kxdx is equivalent to kxdx + 2*pi, so kxdx is
%           on the folded complex plane (ie, the surface of a cylinder).
%           There are 2*ny unique solutions of kxdx, and whenever kxdx is a
%           solution, -kxdx is also a solution. Here, we choose the ny solutions
%           such that
%           1) When k0dx is real, we have
%               Propagating channels: 0 < Re(kxdx) < pi, Im(kxdx) = 0.
%               Evanescent channels: Im(kxdx) >= 0, mod(Re(kxdx),pi) = 0.
%           2) When k0dx is complex, we analytically continue the above choice
%           onto the complex-frequency plane. Specifically, we pick the kxdx
%           that is continuously connected to the one with real k0dx through a
%           vertical line in the complex-(k0dx^2) plane.
%       channels.L.sqrt_mu (1-by-N_prop row vector):
%           Square root of the normalized longitudinal group velocity of the
%           propagating channels, sqrt_mu = sqrt(sin(kxdx)). The longitudinal
%           group velocity is v_g = (sin(kxdx)/k0dx)*(c/epsilon_L).
%       channels.L.ind_prop (1-by-N_prop integer row vector):
%           Indices of the N_prop propagating channels among all ny channels.
%       channels.L.ind_prop_conj (1-by-N_prop integer row vector, optional):
%           A permutation vector that switches one propagating channel with one
%           having a complex-conjugated transverse profile. If
%               kydx_prop = channels.kydx(channels.L.ind_prop),
%           then
%               channels.fun_phi(kydx_prop(channels.L.ind_prop_conj))
%           equals
%               conj(channels.fun_phi(kydx_prop)).
%           For periodic BC, this flips the sign of ky. For Dirichlet and
%           Neumann BC, fun_phi is real, so there is no permutation. If r is the
%           reflection matrix computed using the Fisher-Lee relation,
%           r(channels.L.ind_prop_conj,:) will be symmetric, following from the
%           symmetry of the wave operator (ie, Lorentz reciprocity).
%           For Bloch periodic BC with ka != 0, complex conjugation maps ky to
%           -ky which only exists for a different Bloch BC with -ka, so such
%           permutation does not exist, and ind_prop_conj is not given.

% Check input parameters
if ~(isscalar(ny) && isnumeric(ny) && isreal(ny) && (round(ny)==ny) && ny>0); error('ny must be a positive integer scalar'); end
if ~((ischar(BC) && isvector(BC)) || ((isstring(BC) || isnumeric(BC)) && isscalar(BC)))
    error('Boundary condition BC has to be a character vector or string, or numeric scalar (for Bloch periodic BC)');
end
if ~(isscalar(k0dx) && isnumeric(k0dx)); error('k0dx must be a numeric scalar'); end
if ~(isscalar(epsilon_L) && isnumeric(epsilon_L)); error('epsilon_L must be a numeric scalar'); end
if nargin == 4 || isempty(epsilon_R) || isnan(epsilon_R)
    two_sided = false;
else
    if ~(isscalar(epsilon_R) && isnumeric(epsilon_R)); error('epsilon_R must be a numeric scalar'); end
    two_sided = true;
end

% these are used only for periodic and Bloch periodic boundary conditions; otherwise they stay empty
ka = [];
ind_zero_ky = [];

% handle periodic and Bloch periodic boundary conditions
if strcmpi(BC, 'Bloch')
    error('To use Bloch periodic boundary condition, set BC to k*a where k is the Bloch wave number and a is the periodicity');
elseif isnumeric(BC)
    ka = BC;
    BC = 'Bloch';
    % ka must be real for channels.fun_phi(channels.kydx) to be unitary
    if ~isreal(ka)
        error('k*a = %g + 1i*%g is a complex number; has to be real for a complete orthonormal transverse basis.', real(ka), imag(ka));
    end
elseif strcmpi(BC, 'periodic')
    ka = 0;
    BC = 'Bloch';
end

% f = [f(1), ..., f(ny)].'; 
% For periodic and Bloch periodic BC, we order channels.kydx such that it increases monotonically from negative to positive
% For other BC, ky >= 0, and we order channels.kydx such that it increases monotonically from smallest to largest

% Transverse modes (form a complete basis and are independent of epsilon_L/R)
if strcmpi(BC, 'Bloch')
    % f(ny+1) = f(1)*exp(1i*ka); f(0) = f(ny)*exp(-1i*ka)
    % The transverse mode index where kydx = ka/ny
    if mod(ny,2) == 1
        ind_zero_ky = round((ny+1)/2);
    else
        ind_zero_ky = round(ny/2);
    end
    channels.kydx = (ka/ny) + ((1:ny)-ind_zero_ky)*(2*pi/ny);
    % Normalized transverse mode profile: phi_{m,a} = exp(i*(m-m0)*kydx(a))/sqrt(ny)
    % y=0 is centered at m=m0, and we let m0=(ny+1)/2 in the middle
    channels.fun_phi = @(kydx) exp(((1:ny).')*(1i*kydx))/sqrt(ny);
elseif strcmpi(BC, 'Dirichlet') || strcmpi(BC, 'PEC') % PEC on both sides
    % f(0) = f(ny+1) = 0
    channels.kydx = (1:ny)*(pi/(ny+1));
    % Normalized transverse mode profile: phi_{m,a} = sin(m*kydx(a))*sqrt(2/(ny+1))
    channels.fun_phi = @(kydx) sin(((1:ny).')*kydx)*sqrt(2/(ny+1)); 
elseif strcmpi(BC, 'Neumann') || strcmpi(BC, 'PMC') % PMC on both sides
    % f(0) = f(1), f(ny+1) = f(ny)
    channels.kydx = ((1:ny)-1)*(pi/ny);
    % Normalized transverse mode profile:
    % When kydx == 0: phi_{m,a} = sqrt(1/ny)
    % When kydx != 0: phi_{m,a} = cos((m-0.5)*kydx(a))*sqrt(2/ny)
    % We subtract (~kydx)*(1-sqrt(1/2)) from the cos() which is nonzero only when kydx=0
    channels.fun_phi = @(kydx) (cos(((0.5:ny).')*kydx)-((~kydx)*(1-sqrt(1/2))))*sqrt(2/ny); 
elseif strcmpi(BC, 'DirichletNeumann') || strcmpi(BC, 'PECPMC') % PEC on the low side, PMC on the high side
    % f(0) = 0, f(ny+1) = f(ny)
    channels.kydx = (0.5:ny)*(pi/(ny+0.5));
    % Normalized transverse mode profile: phi_{m,a} = sin(m*kydx(a))*sqrt(2/(ny+0.5))
    channels.fun_phi = @(kydx) sin(((1:ny).')*kydx)*sqrt(2/(ny+0.5)); 
elseif strcmpi(BC, 'NeumannDirichlet') || strcmpi(BC, 'PMCPEC') % PMC on the low side, PEC on the high side
    % f(0) = f(1), f(ny+1) = 0
    channels.kydx = (0.5:ny)*(pi/(ny+0.5));
    % Normalized transverse mode profile: phi_{m,a} = cos((m-0.5)*kydx(a))*sqrt(2/(ny+0.5))
    channels.fun_phi = @(kydx) cos(((0.5:ny).')*kydx)*sqrt(2/(ny+0.5)); 
else
    error('Boundary condition BC = %s is not an eligible option', BC);
end

% Longitudinal properties for homogeneous space on the left (kxdx, sqrt_mu, number of propagating channels, etc; depends on epsilon_L/R)
channels.L = setup_longitudinal((k0dx^2)*epsilon_L, channels.kydx, ka, ind_zero_ky);

% Homogeneous space on the right
if two_sided
    if epsilon_R == epsilon_L
        channels.R = channels.L;
    else
        channels.R = setup_longitudinal((k0dx^2)*epsilon_R, channels.kydx, ka, ind_zero_ky);
    end
elseif nargin == 5
    % create channels.R.N_prop as this variable may be accessed even in one-sided systems
    channels.R.N_prop = NaN;
end

end


function side = setup_longitudinal(k0dx2_epsilon, kydx, ka, ind_zero_ky)
% Returns a structure 'side'. See comments at the beginning of this file for more info.

% finite-difference dispersion for homogeneous space: k0dx2_epsilon = 4*sin^2(kxdx/2) + 4*sin^2(kydx/2)

% sin_kxdx_over_two_sq = sin^2(kxdx/2)
sin_kxdx_over_two_sq = 0.25*k0dx2_epsilon - sin(kydx/2).^2;

% Normalized longitudinal wave number
% asin(sqrt(z)) has two branch points (at z=0 and z=1) and with the branch cuts going outward on the real-z axis; we will address the branch choice below
% Note kxdx is only defined up to modulo 2*pi (ie, kxdx is equivalent to kxdx + 2*pi, kxdx - 2*pi, etc) because sin(kxdx) and exp(1i*kxdx) are both invariant under 2pi shifts.
side.kxdx = 2*asin(sqrt(sin_kxdx_over_two_sq));

% Indices of the propagating channels
% When k0dx2_epsilon is real, these are indicies of the channels with real-valued kxdx
% When k0dx2_epsilon is complex, these are indicies of the channels we consider "propagating-like"; they have complex kxdx with 0 < real(kxdx) < pi. When k0dx2_epsilon is tuned to a real number continuously, this set continuously becomes that at real k0dx2_epsilon.
side.ind_prop = find((real(sin_kxdx_over_two_sq) > 0) & (real(sin_kxdx_over_two_sq) < 1));

% Number of propagating channels
side.N_prop = length(side.ind_prop);
if side.N_prop==0 && length(kydx)==1; side.ind_prop = zeros(1,0); end  % a rare possibility, but in this case ind_prop would be zeros(0,0) while it should be zeros(1,0)

% Here we address the sign choice of kxdx, namely its branch
% When k0dx2_epsilon is real, we choose the sign of kxdx such that:
% 1) 0 < kxdx < pi for propagating channels (where kxdx is real)
% 2) Im(kxdx) >= 0 for evanescent channels
% Using the correct sign is important when we build the retarded Green's function of the semi-infinite homogeneous space.
% The default branch choice of asin(sqrt(z)) returns the correct sign for the most part, except when z > 1. We need to flip the sign of those (which can only occur if k0dx2_epsilon > 4).
% When k0dx2_epsilon is complex-valued, it is not always possible to unambiguously choose the sign that is "physical", because kxdx will be complex-valued, and the sign we "want" for real(kxdx) and the sign we want for imag(kxdx) may be incompatible.
% What we do with complex k0dx2_epsilon is that we choose the sign for the (complex-valued) kxdx such that when k0dx2_epsilon is tuned to a real number continuously by fixing Re(k0dx2_epsilon) and varying Im(k0dx2_epsilon), the kxdx we choose continuously becomes the "correct" one at real k0dx2_epsilon without crossing any branch cut. To do so, we rotate the two branch cuts of asin(sqrt(z)) by 90 degrees to the lower part of the complex-z plane (ie, the lower part of the complex-k0dx2_epsilon plane), and we pick the branch with the correct sign when k0dx2_epsilon is real. This is implemented by flipping the sign of kxdx for the ones that require flipping.
% Note that we will get a discontinuity whenever k0dx2_epsilon crosses one of those vertical-pointing branch cuts. That is unavoidable.
% The following few lines implement the "flipping".
if ~isreal(k0dx2_epsilon) || (isreal(k0dx2_epsilon) && k0dx2_epsilon > 4)
    % Note that when imag(sin_kxdx_over_two_sq)=0, flipping is needed for sin_kxdx_over_two_sq>1 but not needed for sin_kxdx_over_two_sq<0.
    ind_flip = find((real(sin_kxdx_over_two_sq)<0 & imag(sin_kxdx_over_two_sq)<0) | (real(sin_kxdx_over_two_sq)>1 & imag(sin_kxdx_over_two_sq)<=0));
    side.kxdx(ind_flip) = -side.kxdx(ind_flip);
end

% Square root of the normalized longitudinal group velocity, sqrt(sin(kxdx)), for the propagating channels
% When k0dx2_epsilon is real, sqrt_mu is also real. When k0dx2_epsilon is complex, sqrt_mu is also complex.
side.sqrt_mu = sqrt(sin(side.kxdx(side.ind_prop)));

% Permutation that switches one propagating channel with one having a complex-conjugated transverse profile.
if isempty(ka)
    % For Dirichlet and Neumann BC, fun_phi is real, so no permutation needed
    side.ind_prop_conj = 1:side.N_prop;
else
    if ka == 0
        % For periodic boundary condition, complex conjugation switches ky and -ky
        if ismember(ind_zero_ky, side.ind_prop) || (mod(side.N_prop,2)==0)
            % Simply flip the ordering
            side.ind_prop_conj = side.N_prop:-1:1;
        else
            % The last channel has -ky equal to ky due to aliasing so should not be flipped
            side.ind_prop_conj = [(side.N_prop-1):-1:1, side.N_prop];
        end
    end
end

end
