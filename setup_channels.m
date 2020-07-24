function channels = setup_channels(ny, BC, kdx, epsilon_L, epsilon_R)
%SETUP_CHANNELS    Set up channels of the homogeneous spaces.
%   Returns a structure containing properties of the propagating and evanescent
%   channels in the homogeneous spaces on the left and right.
%
%   === Input parameters ===
%   ny (positive integer scalar):
%       Number of grid points in the transverse (y) direction.
%   BC (character vector or string, case insensitive):
%       Boundary condition in the transverse direction y; can be 'periodic' or
%       'Dirichlet'. Periodic BC means E_z(n,m+ny)=E_z(n,m). Dirichlet BC means
%       E_z(n,m=0)=E_z(n,m=ny+1)=0.
%   kdx (numeric scalar, real or complex):
%       Normalized frequency k*dx.
%   epsilon_L (numeric scalar, real or complex):
%       Relative permittivity of the homogeneous space on the left.
%   epsilon_R (numeric scalar, real or complex):
%       Relative permittivity of the homogeneous space on the right.
%
%   === Returns ===
%   channels (structure):
%       channels.L and channels.R (structue):
%           Structures containing properties specific to the left (L) and right
%           (R) sides.
%       channels.L.N_prop (integer scalar):
%           Number of propagating channels.
%       channels.L.kxdx (1-by-ny complex row vector):
%           Normalized longitudinal wave number kx*dx for all channels. Due to
%           the discretization, kxdx is equivalent to kxdx + 2*pi, so kxdx is
%           on the folded complex plane (ie, the surface of a cylinder).
%           Thereare 2*ny unique solutions of kxdx, and whenever kxdx is a
%           solution, -kxdx is also a solution. Here, we choose the ny solutions
%           such that
%           1) When kdx is real, we have
%               Propagating channels: 0 <= Re(kxdx) <= pi, Im(kxdx) = 0.
%               Evanescent channels: Im(kxdx) > 0, mod(Re(kxdx),pi) = 0.
%           2) When kdx is complex, we analytically continue the above choice
%           onto the complex-frequency plane. Specifically, we pick the kxdx
%           that is continuously connected to the one with real kdx through a
%           vertical line in the complex-(kdx^2) plane.
%       channels.L.sqrt_vg (1-by-N_prop row vector):
%           Square root of the normalized longitudinal group velocity of the
%           propagating channels, sqrt_vg = sqrt(sin(kxdx)). The longitudinal
%           group velocity is v_g = (sin(kxdx)/kdx)*(c/epsilon_L).
%       channels.L.ind_prop (1-by-N_prop integer row vector):
%           Indices of the N_prop propagating channels among all ny channels.
%       channels.L.ind_prop_conj (1-by-N_prop integer row vector):
%           A permutation vector that switches one propagating channel with one
%           having a complex-conjugated transverse profile. For periodic BC,
%           this flips the sign of the propagation angle. For Dirichlet BC,
%           there is no permutation. This permutation can be used to ensure
%           symmetry of the scattering matrix required by Lorentz reciprocity.
%           For example, if r is the reflection matrix, then r itself is not
%           necessarily symmetric, but r(channels.L.ind_prop_conj,:) should be.
%       channels.kydx (1-by-ny real row vector):
%           Normalized transverse wave number ky*dx for all ny channels. They
%           are real-valued and are ordered from small to large.
%       channels.fun_chi (function_handle):
%           A function that, given one element of kydx as the input, returns its
%           normalized transverse field profile as an ny-by-1 column vector. The
%           transverse modes form a complete and orthonormal set, so the
%           ny-by-ny matrix chi = fun_chi(kydx) is unitary.
%       channels.ind_conj (1-by-ny integer row vector):
%           A permutation vector that switches one channel with one having a
%           complex-conjugated field profile. For periodic BC, this flips the
%           sign of kydx. For Dirichlet BC, there is no permutation.

% Check input parameters
if ~(isscalar(ny) && isnumeric(ny) && isreal(ny) && (round(ny)==ny) && ny>0); error('ny must be a positive integer scalar'); end
if ~((ischar(BC) || isstring(BC)) && isrow(BC)); error('BC must be ''periodic'' or ''Dirichlet'''); end
if ~(isscalar(kdx) && isnumeric(kdx)); error('kdx must be a numeric scalar'); end
if ~(isscalar(epsilon_L) && isnumeric(epsilon_L)); error('epsilon_L must be a numeric scalar'); end
if ~(isscalar(epsilon_R) && isnumeric(epsilon_R)); error('epsilon_R must be a numeric scalar'); end

% Transverse modes (form a complete basis and are independent of epsilon_L/R)
if strcmpi(BC, 'Dirichlet')
    ind_zero_ky = 0; % dummy variable; won't be used
    channels.ind_conj = 1:ny; % No permutation needed
    % Order the transverse modes such that ky increases monotonically from smallest to largest
    channels.kydx = (1:ny)*(pi/(ny+1));
    % Normalized transverse mode profile: chi_{m,a} = sin(m*kydx(a))*sqrt(2/(ny+1))
    channels.fun_chi = @(kydx) sin(((1:ny).')*kydx)*sqrt(2/(ny+1)); 
else % periodic
    if mod(ny,2) == 1
        % The transverse mode index where kydx = 0
        ind_zero_ky = round((ny+1)/2);
        % Permutation: simply flip the ordering
        channels.ind_conj = ny:-1:1;
    else
        % The transverse mode index where kydx = 0
        ind_zero_ky = round(ny/2);
        % Permutation: the last mode has -ky equal to ky due to aliasing so should not be flipped
        channels.ind_conj = [(ny-1):-1:1, ny];
    end
    % Order the transverse modes such that ky increases monotonically from negative to positive
    channels.kydx = ((1:ny)-ind_zero_ky)*(2*pi/ny);
    % Normalized transverse mode profile: chi_{m,a} = exp(i*m*kydx(a))/sqrt(ny)
    channels.fun_chi = @(kydx) exp(((1:ny).')*(1i*kydx))/sqrt(ny);
end

% Longitudinal properties for homogeneous space on the left (kxdx, sqrt_vg, number of propagating channels, etc; depends on epsilon_L/R)
channels.L = setup_longitudinal(BC, (kdx^2)*epsilon_L, 'left', channels.kydx, ind_zero_ky);

% Check that the permutation vectors are consistent
if ~isequal(channels.ind_conj(channels.L.ind_prop(channels.L.ind_prop_conj)), channels.L.ind_prop); error('channel ordering incorrect'); end

% Homogeneous space on the right
if epsilon_R == epsilon_L
    channels.R = channels.L;
else
    channels.R = setup_longitudinal(BC, (kdx^2)*epsilon_R, 'right', channels.kydx, ind_zero_ky);
    if ~isequal(channels.ind_conj(channels.R.ind_prop(channels.R.ind_prop_conj)), channels.R.ind_prop); error('channel ordering incorrect'); end
end

end


function side = setup_longitudinal(BC, kdx2_epsilon, side_str, kydx, ind_zero_ky)
% Returns a structure 'side'. See comments at the beginning of this file for more info.

% cos(kxdx) from the disperision relation for homogeneous space in the finite-difference wave equation
cos_kxdx = 2 - 0.5*kdx2_epsilon - cos(kydx);

%if min(abs(cos_kxdx-1)) < 1e-12 || min(abs(cos_kxdx+1)) < 1e-12
%    warning('%s channels: exist channel(s) at the cutoff between propagating and evanescent channels', side_str);
%end

% Indices of the propagating channels
% When kdx2_epsilon is real, these are indicies of the channels with real-valued kxdx
% When kdx2_epsilon is complex, these are indicies of the channels we consider "propagating-like"; they have complex kxdx with 0 < real(kxdx) < pi. When kdx2_epsilon is tuned to a real number continuously, this set continuously becomes that at real kdx2_epsilon.
side.ind_prop = find(abs(real(cos_kxdx)) <= 1);
side.N_prop = length(side.ind_prop);
if side.N_prop==0 && length(cos_kxdx)==1; side.ind_prop = zeros(1,0); end  % a rare possibility, but in this case ind_prop would be zeros(0,0) while it should be zeros(1,0)

% Normalized longitudinal wave number
% Note that acos() in MATLAB returns the principal branch of acos(z), which has two branch points (at z=1 and z=-1) and with the branch cuts going outward on the real axis (at |z|>1).
side.kxdx = acos(cos_kxdx);

% Note kxdx is only defined up to modulo 2*pi (ie, kxdx is equivalent to kxdx + 2*pi, kxdx - 2*pi, etc), because sin(kxdx) and exp(1i*kxdx) are both invariant under 2pi shifts.
% However, exp(1i*kxdx) and exp(1i*(-kxdx)) are very different. So, the sign choice of kxdx is important.
% When kdx2_epsilon is real, we want to choose the sign of kxdx such that:
% 1) 0 <= kxdx <= pi for propagating channels (where kxdx is real)
% 2) Im(kxdx) > 0 for evanescent channels
% The principal branch of acos returns the correct sign for the most part, except when cos_kxdx < -1, for which acos() returns the branch with Im(kxdx) < 0, which is not what we want. We need to flip the sign of those (which can only occur if kdx2_epsilon > 4).
% When kdx2_epsilon is complex-valued, it is not always possible to unambiguously choose the sign that is "physical", because kxdx will be complex-valued, and the sign we "want" for real(kxdx) and the sign we want for imag(kxdx) may be incompatible.
% What we do with complex kdx2_epsilon is that we choose the sign for the (complex-valued) kxdx such that when kdx2_epsilon is tuned to a real number continuously, the kxdx we choose continuously becomes the "correct" one at real kdx2_epsilon. To do so, we rotate the two branch cuts of acos() by 90 degrees to the upper part of the complex-cos_kxdx plane (ie, the lower part of the complex-kdx2_epsilon plane) so that no branch-crossing occurs when kdx2_epsilon approaches the real axis, and we pick the branch that contains the "correct" sign when kdx2_epsilon is real. This is implemented by flipping the sign of kxdx, for the ones that require flipping.
% Note that we will get a discontinuity whenever kdx2_epsilon crosses one of those vertical-pointing branch cuts. That is unavoidable.
% The following few lines implement the "flipping".
if ~isreal(kdx2_epsilon) || (isreal(kdx2_epsilon) && kdx2_epsilon > 4)
    % indicies where flipping is necessary to go from the principal branch of acos to the desired branch described above.
    % Note that when imag(cos_kxdx)=0, flipping is needed for cos_kxdx<-1 but not needed for cos_kxdx>1, as described above.
    ind_flip = find((real(cos_kxdx)>1 & imag(cos_kxdx)>0) | (real(cos_kxdx)<-1 & imag(cos_kxdx)>=0));
    side.kxdx(ind_flip) = -side.kxdx(ind_flip);
end

% Square root of the normalized longitudinal group velocity, sqrt(sin(kxdx)), for the propagating channels
% When kdx2_epsilon is real, sqrt_vg is also real. When kdx2_epsilon is complex, sqrt_vg is also complex.
side.sqrt_vg = sqrt(sin(side.kxdx(side.ind_prop)));

% Permutation that switches one propagating channel with one having a complex-conjugated transverse profile.
switch lower(BC)
    case 'periodic'
        if ismember(ind_zero_ky, side.ind_prop) || (mod(side.N_prop,2)==0)
            % Simply flip the ordering
            side.ind_prop_conj = side.N_prop:-1:1;
        else
            % The last channel has -ky equal to ky due to aliasing so should not be flipped
            side.ind_prop_conj = [(side.N_prop-1):-1:1, side.N_prop];
        end
    case 'dirichlet'
        % No permutation needed
        side.ind_prop_conj = 1:side.N_prop;
    otherwise
        error('boundary condition ''%s'' is not supported', string(BC));
end

end
