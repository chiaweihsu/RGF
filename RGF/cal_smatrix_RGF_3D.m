function [smatrix, channels] = cal_smatrix_RGF_3D(syst, out, in)
%CAL_SMATRIX_RGF_3D   Scattering matrix from the recursive Green's function method.
%   [smatrix, channels] = CAL_SMATRIX_RGF_3D(syst, out, in) returns portions of
%   the scattering matrix with output and input channels specified by out and
%   in, for the system specified by syst. The system is the scalar wave equation
%   in 3D:
%       [(d/dx)^2 + (d/dy)^2 + (d/dz)^2 + k^2*epsilon(x,y,z)]*psi(x,y,z) = 0,
%   discretized using center difference on a cubic grid with grid size dx, with
%   homogeneous spaces on the left and on the right.
%
%   === Input parameters ===
%   syst (structure):
%       syst.epsilon (real or complex matrix):
%           Relative permittivity of the scattering region, with size [nx,ny,nz]
%            = size(syst.epsilon). Absorption can be included through a positive
%           imaginary part. Linear gain can be included with a negative
%           imaginary part.
%       syst.epsilon_L (numeric scalar, real or complex):
%           Relative permittivity of the homogeneous space on the left.
%       syst.epsilon_R (numeric scalar, real or complex):
%           Relative permittivity of the homogeneous space on the right.
%       syst.k0dx (numeric scalar, real or complex):
%           Normalized frequency k0*dx = (2*pi/lambda)*dx.
%       syst.xyBC (character vector or string (case insensitive) or real vector):
%           Boundary condition in the x and y directions.
%           For character inputs, the options are:
%             'periodic': f(n+nx,m) = f(n,m+ny) = f(n,m)
%             'Dirichlet' or 'PEC': f(0,m) = f(nx+1,m) = f(n,0) = f(n,ny+1) = 0
%           When xyBC is a real vector, the Bloch periodic boundary condition is used
%           with f(n+nx,m) = f(n,m)*exp(1i*kxa), f(n,m+ny) = f(n,m)*exp(1i*kya).
%           Here, xyBC = [kxa, kya]
%   out (cell array or structure or character vector):
%       The output channels of interest. Possible choices:
%       A cell array of character vectors: this cell array may contain 'left'
%           or 'L' to specify all output channels on the left, and/or 'right' or
%           'R' to specify all output channels on the right.
%       A structure: this structure may contain fields 'ind_out_L' and/or
%           'ind_out_R', which are vectors whose elements are the indices of
%           outgoing channels on the left and/or right.
%   in (cell array or structure):
%       The input channels of interest. Possible choices:
%       A cell array of character vectors: this cell array may contain 'left'
%           or 'L' to specify all input channels on the left, and/or 'right' or
%           'R' to specify all input channels on the right.
%       A structure: this structure may contain fields 'ind_in_L' and/or
%           'ind_in_R', which are vectors whose elements are the indices of
%           incoming channels on the left and/or right.
%
%   === Returns ===
%   smatrix (matrix):
%       A block of the scattering matrix with output channels (rows) and input
%       channels (columns) specified by 'out' and 'in'.
%       == Ordering ==
%       Channels on the left are always ordered before channels on the right.
%       When 'out' or 'in' is a cell array, the channel indices are
%       1:channels.L.N_prop and 1:channels.R.N_prop for each side. When 'out' or
%       'in' is a structure, the channel indicies are those given by in.ind_in_L,
%       etc.
%       == Channels ==
%       Propagating channel a on the left has wave vector [channels.kxdx(ind),
%       channels.kydx(ind), channels.L.kzdx(ind)] and transverse profile
%       channels.fun_phi(channels.kxdx(ind), channels.kydx(ind)), where
%       ind=channels.L.ind_prop(a).
%       == Reference planes ==
%       The reference planes for the scattering matrix elements are at one pixel
%       outside the scattering region, at n=0 and n=nx+1.
%   channels (structure):
%       A structure returned by function setup_channels_3D(). It contains
%       properties of channels of the homogeneous spaces on the left and right.
%       Type "help setup_channels_3D" for more information.
%
%   See also: setup_channels_3D

%% Check validity of the input arguments

% Check input parameters
if ~(isstruct(syst) && numel(syst)==1); error('syst must be a structure'); end
if ~isfield(syst, 'epsilon');   error('syst.epsilon must be provided'); end
if ~isfield(syst, 'epsilon_L'); error('syst.epsilon_L must be provided'); end
if ~isfield(syst, 'epsilon_R'); error('syst.epsilon_R must be provided'); end
if ~isfield(syst, 'k0dx');      error('syst.k0dx must be provided'); end
if ~isfield(syst, 'xyBC');        error('syst.xyBC must be provided'); end
if ~((ndims(syst.epsilon)==3 || ndims(syst.epsilon)==2) && isnumeric(syst.epsilon)); error('syst.epsilon must be a 3D array'); end
if ~(isscalar(syst.epsilon_L) && isnumeric(syst.epsilon_L)); error('syst.epsilon_L must be a numeric scalar'); end
if ~(isscalar(syst.epsilon_R) && isnumeric(syst.epsilon_R)); error('syst.epsilon_R must be a numeric scalar'); end
if ~(isscalar(syst.k0dx) && isnumeric(syst.k0dx));           error('syst.k0dx must be a numeric scalar'); end

% Number of grid points
[nx, ny, nz] = size(syst.epsilon);
nxy = nx*ny;
if nxy==0; error('need at least one site along the transverse directions (x,y)'); end

%% Build input matrix B and output matrix C (some computations, but not much)

% Set up the homogeneous-space channels on the two sides
channels = setup_channels_3D(nx, ny, syst.xyBC, syst.k0dx, syst.epsilon_L, syst.epsilon_R);
N_prop_L = channels.L.N_prop; N_prop_R = channels.R.N_prop;

%fprintf('nx = %d; ny = %d; nz = %d; N_prop = %d, %d\n', nx, ny, nz, N_prop_L, N_prop_R);

% Input channels, in row vectors
ind_in_L = zeros(1, 0);
ind_in_R = zeros(1, 0);
if isstruct(in)
    % Take the user-specified lists of input channels from left and/or from right
    if isfield(in, 'ind_in_L')
        ind_in_L = reshape(in.ind_in_L, 1, []);
        if ~isempty(ind_in_L) && ~(isnumeric(ind_in_L) && isequal(round(ind_in_L),ind_in_L) && min(ind_in_L)>0 && max(ind_in_L)<=N_prop_L)
            error('in.ind_in_L, when specified, must be an array of positive integers not exceeding N_prop_L=%d', N_prop_L);
        end
    end
    if isfield(in, 'ind_in_R')
        ind_in_R = reshape(in.ind_in_R, 1, []);
        if ~isempty(ind_in_R) && ~(isnumeric(ind_in_R) && isequal(round(ind_in_R),ind_in_R) && min(ind_in_R)>0 && max(ind_in_R)<=N_prop_R)
            error('in.ind_in_R, when specified, must be an array of positive integers not exceeding N_prop_R=%d', N_prop_R);
        end
    end
elseif iscell(in)
    % Take all input channels on the specified side(s)
    for k = 1:numel(in)
        element = in{k};
        if ~((ischar(element) || isstring(element)) && isrow(element))
            warning('element in{%d} has wrong format and is ignored; should be one of {''left'', ''right'', ''L'', ''R''}', k);
        end
        switch char(element)
            case {'left','L'}
                ind_in_L = 1:N_prop_L;
            case {'right','R'}
                ind_in_R = 1:N_prop_R;
            otherwise
                warning('element in{%d} = %s ignored; should be one of {''left'', ''right'', ''L'', ''R''}', k, char(element));
        end
    end
else
    error('input argument ''in'' must be a structure or a cell array')
end

% Output channels, in row vectors
ind_out_L = zeros(1, 0);
ind_out_R = zeros(1, 0);
if isstruct(out)
    % Take the user-specified lists of output channels from left and/or from right
    if isfield(out, 'ind_out_L')
        ind_out_L = reshape(out.ind_out_L, 1, []);
        if ~isempty(ind_out_L) && ~(isnumeric(ind_out_L) && isequal(round(ind_out_L),ind_out_L) && min(ind_out_L)>0 && max(ind_out_L)<=N_prop_L)
            error('out.ind_out_L, when specified, must be an array of positive integers not exceeding N_prop_L=%d', N_prop_L);
        end
    end
    if isfield(out, 'ind_out_R')
        ind_out_R = reshape(out.ind_out_R, 1, []);
        if ~isempty(ind_out_R) && ~(isnumeric(ind_out_R) && isequal(round(ind_out_R),ind_out_R) && min(ind_out_R)>0 && max(ind_out_R)<=N_prop_R)
            error('out.ind_out_R, when specified, must be an array of positive integers not exceeding N_prop_R=%d', N_prop_R);
        end
    end
elseif iscell(out)
    % Take all output channels on the specified side(s)
    for k = 1:numel(out)
        element = out{k};
        if ~((ischar(element) || isstring(element)) && isrow(element))
            warning('element out{%d} has wrong format and is ignored; should be one of {''left'', ''right'', ''L'', ''R''}', k);
        end
        switch char(element)
            case {'left','L'}
                ind_out_L = 1:N_prop_L;
            case {'right','R'}
                ind_out_R = 1:N_prop_R;
            otherwise
                warning('element out{%d} = %s ignored; should be one of {''left'', ''right'', ''L'', ''R''}', k, char(element));
        end
    end
else
    error('input argument ''out'' must be a structure or a cell array')
end

N_in_L = numel(ind_in_L); N_in_R = numel(ind_in_R); N_in_tot = N_in_L + N_in_R;
N_out_L = numel(ind_out_L); N_out_R = numel(ind_out_R); N_out_tot = N_out_L + N_out_R;

% No computation needed if the requested S-matrix is empty
if (N_in_tot==0 || N_out_tot==0)
    smatrix = zeros(N_out_tot, N_in_tot);
    return
end

% Build nxy-by-nxy unitary matrix phi where the a-th column is the a-th transverse mode reshaped into a column vector
phi = zeros(nxy, nxy);
for a = 1:nxy
    phi(:,a) = reshape(channels.fun_phi(channels.kxdx(a), channels.kydx(a)), [], 1);
end

% Build input matrix B and output matrix C; the -2i prefactor will be multiplied at the end
% B and C at one pixel outside the scattering region (n=0 and n=nx+1), as dense matrices
B_L = phi(:,channels.L.ind_prop(ind_in_L))*spdiags(reshape(channels.L.sqrt_mu(ind_in_L),[],1), 0, N_in_L, N_in_L);
B_R = phi(:,channels.R.ind_prop(ind_in_R))*spdiags(reshape(channels.R.sqrt_mu(ind_in_R),[],1), 0, N_in_R, N_in_R);
% Note that when the frequency k0dx is complex, sqrt_mu is also complex, and we should take the conjugate transpose of phi but not of sqrt_mu
C_L = spdiags(reshape(channels.L.sqrt_mu(ind_out_L),[],1), 0, N_out_L, N_out_L)*(phi(:,channels.L.ind_prop(ind_out_L))');
C_R = spdiags(reshape(channels.R.sqrt_mu(ind_out_R),[],1), 0, N_out_R, N_out_R)*(phi(:,channels.R.ind_prop(ind_out_R))');

% Retarded Green's function G0 of a semi-infinite homogeneous space, evaluated at the surface (just before the space is terminated)
G0_L = phi*spdiags(exp(1i*channels.L.kzdx(:)), 0, nxy, nxy)*(phi');
if syst.epsilon_R == syst.epsilon_L
    G0_R = G0_L;
else
    G0_R = phi*spdiags(exp(1i*channels.R.kzdx(:)), 0, nxy, nxy)*(phi');
end

clear phi

% Bloch periodic BC
if isnumeric(syst.xyBC)
    xBC = syst.xyBC(1);
    yBC = syst.xyBC(2);
else
    xBC = syst.xyBC;
    yBC = syst.xyBC;
end

% second derivatives in x and y directions
[laplacian_x] = build_laplacian_1d(nx, xBC);
[laplacian_y] = build_laplacian_1d(ny, yBC);

% Finite-difference wave operator for one slice, without the (k0dx)^2*epsilon term
% use Kronecker outer product to go from 1D to 2D
% adds 2 for the laplacian in z direction
A0 = -kron(speye(ny), laplacian_x) - kron(laplacian_y, speye(nx)) + 2*speye(nxy);

k0dx2 = (syst.k0dx)^2;

%% Main computation

% Iterate through by adding slices
% To avoid unnecessary steps, we pick the scheme based on which parts of the S-matrix are requested
if N_in_R==0 && N_out_R==0
    % input and output both on the left; loop from right to left
    G_LL = G0_R;
    for n = nz:-1:1
        A_nn = A0 - spdiags(k0dx2*reshape(syst.epsilon(:,:,n),[],1), 0, nxy, nxy);
        G_LL = inv(A_nn - G_LL);
    end
    G_LL = (eye(nxy) - G0_L * G_LL) \ G0_L;
    smatrix = (-2i)*(C_L * (G_LL * B_L));
elseif N_in_R==0
    % input from the left; loop from right to left
    G_LL = G0_R; G_RL = G0_R;
    for n = nz:-1:1
        A_nn = A0 - spdiags(k0dx2*reshape(syst.epsilon(:,:,n),[],1), 0, nxy, nxy);
        G_LL = inv(A_nn - G_LL);
        G_RL = - G_RL * G_LL;
    end
    G_LL = (eye(nxy) - G0_L * G_LL) \ G0_L;
    G_RL = - G_RL * G_LL;
    smatrix = (-2i)*[(C_L * G_LL) * B_L; C_R * (G_RL * B_L)]; % also works when N_out_L==0
elseif N_in_L==0 && N_out_L==0
    % input and output both on the right; loop from left to right
    G_RR = G0_L;
    for n = 1:nz
        A_nn = A0 - spdiags(k0dx2*reshape(syst.epsilon(:,:,n),[],1), 0, nxy, nxy);
        G_RR = inv(A_nn - G_RR);
    end
    G_RR = (eye(nxy) - G0_R * G_RR) \ G0_R;
    smatrix = (-2i)*(C_R * (G_RR * B_R));
elseif N_in_L==0
    % input from the right; loop from left to right
    G_RR = G0_L; G_LR = G0_L;
    for n = 1:nz
        A_nn = A0 - spdiags(k0dx2*reshape(syst.epsilon(:,:,n),[],1), 0, nxy, nxy);
        G_RR = inv(A_nn - G_RR);
        G_LR = - G_LR * G_RR;
    end
    G_RR = (eye(nxy) - G0_R * G_RR) \ G0_R;
    G_LR = - G_LR * G_RR;
    smatrix = (-2i)*[C_L * (G_LR * B_R); (C_R * G_RR) * B_R]; % also works when N_out_R==0
else
    % general case: input from both sides; loop from left to right
    % Initial step: retarded Green's function for an semi-infinite homogeneous space on the left
    G_RR = G0_L; G_RL = G0_L; G_LL = G0_L; G_LR = G0_L;

    % Main loop of RGF: attach the scattering region slice by slice
    for n = 1:nz
        A_nn = A0 - spdiags(k0dx2*reshape(syst.epsilon(:,:,n),[],1), 0, nxy, nxy);
        G_RR = inv(A_nn - G_RR);
        G_RL = - G_RR * G_RL;
        G_LL =   G_LL - G_LR * G_RL;
        G_LR = - G_LR * G_RR;
    end

    % Final step: attach homogeneous space on the right
    G_RR = (eye(nxy) - G0_R * G_RR) \ G0_R;
    G_RL = - G_RR * G_RL;
    G_LL =   G_LL - G_LR * G_RL;
    G_LR = - G_LR * G_RR;

    % Fisher-Lee relation for S = [[r;t],[tp,rp]]; the Kronecker delta term will be included later
    % Multiply to the right first, since typically the number of input channels equals or is less than the number of output channels
    smatrix = (-2i)*[[C_L * (G_LL * B_L); C_R * (G_RL * B_L)], [C_L * (G_LR * B_R); C_R * (G_RR * B_R)]];
end
if issparse(smatrix); smatrix = full(smatrix); end

%% Wrapping up

% Subtract the Kronecker delta term in the Fisher-Lee relation; take the input/output channels with the same indices
[~, ind_out_L_temp, ind_in_L_temp] = intersect(ind_out_L, ind_in_L); % note: the temp indices are column vectors
[~, ind_out_R_temp, ind_in_R_temp] = intersect(ind_out_R, ind_in_R); % note: the temp indices are column vectors
D = sparse([ind_out_L_temp; N_out_L + ind_out_R_temp], [ind_in_L_temp; N_in_L + ind_in_R_temp], 1, N_out_tot, N_in_tot);
smatrix = smatrix - D;

end
