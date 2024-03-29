function [smatrix, channels] = cal_smatrix_RGF(syst, out, in)
%CAL_SMATRIX_RGF   Scattering matrix from the recursive Green's function method.
%   [smatrix, channels] = CAL_SMATRIX_RGF(syst, out, in) returns portions of the
%   scattering matrix with output and input channels specified by out and in,
%   for the system specified by syst. The system is Maxwell's equations for 2D
%   TM modes:
%       [(d/dx)^2 + (d/dy)^2 + k^2*epsilon(x,y)]*E_z(x,y) = 0,
%   discretized using center difference on a square grid with grid size dx, with
%   homogeneous spaces on the left and on the right.
%
%   === Input parameters ===
%   syst (structure):
%       syst.epsilon (real or complex matrix):
%           Relative permittivity of the scattering region, with size [ny,nx] =
%           size(syst.epsilon) where ny>=1, nx>=0. Absorption can be included
%           through a positive imaginary part. Linear gain can be included with
%           a negative imaginary part.
%       syst.epsilon_L (numeric scalar, real or complex):
%           Relative permittivity of the homogeneous space on the left.
%       syst.epsilon_R (numeric scalar, real or complex):
%           Relative permittivity of the homogeneous space on the right.
%       syst.k0dx (numeric scalar, real or complex):
%           Normalized frequency k0*dx = (2*pi/lambda)*dx.
%       syst.yBC (character vector or string (case insensitive) or scalar number):
%           Boundary condition in the y direction.
%           For character inputs, the options are:
%             'periodic': E_z(m+ny,n) = E_z(m,n)
%             'Dirichlet' or 'PEC': E_z(0,n) = E_z(ny+1,n) = 0
%             'Neumann' or 'PMC': E_z(0,n) = E_z(1,n), E_z(ny+1,n) = E_z(ny,n)
%             'DirichletNeumann' or 'PECPMC': E_z(0,n) = 0, E_z(ny+1,n) = E_z(ny,n)
%             'NeumannDirichlet' or 'PMCPEC': E_z(0,n) = E_z(1,n), E_z(ny+1,n) = 0
%           When syst.yBC is a scalar number, the Bloch periodic boundary
%           condition is used with E_z(m+ny,n) = E_z(m,n)*exp(1i*ky*a) where ky
%           is the Bloch wave number. Here, syst.yBC = ky*a = ky*dx*ny.
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
%       Propagating channel a on the left has wavenumbers in x and y being
%       channels.L.kxdx(ind) and channels.kydx(ind) and transverse profile being
%       channels.fun_phi(channels.kydx(ind)), where ind=channels.L.ind_prop(a).
%       == Reference planes ==
%       The reference planes for the scattering matrix elements are at one pixel
%       outside the scattering region, at n=0 and n=nx+1.
%   channels (structure):
%       A structure returned by function setup_channels(). It contains
%       properties of channels of the homogeneous spaces on the left and right.
%       Type "help setup_channels" for more information.
%
%   See also: setup_channels 

%% Check validity of the input arguments

% Check input parameters
if ~(isstruct(syst) && numel(syst)==1); error('syst must be a structure'); end
if ~isfield(syst, 'epsilon');   error('syst.epsilon must be provided'); end
if ~isfield(syst, 'epsilon_L'); error('syst.epsilon_L must be provided'); end
if ~isfield(syst, 'epsilon_R'); error('syst.epsilon_R must be provided'); end
if ~isfield(syst, 'k0dx');      error('syst.k0dx must be provided'); end
if ~isfield(syst, 'yBC');       error('syst.yBC must be provided'); end
if ~(ismatrix(syst.epsilon) && isnumeric(syst.epsilon));     error('syst.epsilon must be a numeric matrix'); end
if ~(isscalar(syst.epsilon_L) && isnumeric(syst.epsilon_L)); error('syst.epsilon_L must be a numeric scalar'); end
if ~(isscalar(syst.epsilon_R) && isnumeric(syst.epsilon_R)); error('syst.epsilon_R must be a numeric scalar'); end
if ~(isscalar(syst.k0dx) && isnumeric(syst.k0dx));           error('syst.k0dx must be a numeric scalar'); end
if ~((ischar(syst.yBC) || isstring(syst.yBC)) && isrow(syst.yBC)); error('syst.yBC must be ''periodic'' or ''Dirichlet'''); end

% Number of grid points in y and x
[ny, nx] = size(syst.epsilon);
if ny==0; error('need at least one site along the transverse (y) direction'); end

%% Build input matrix B and output matrix C (some computations, but not much)

% Set up the homogeneous-space channels on the two sides
channels = setup_channels(ny, syst.yBC, syst.k0dx, syst.epsilon_L, syst.epsilon_R);
N_prop_L = channels.L.N_prop; N_prop_R = channels.R.N_prop;

%fprintf('nx = %d; ny = %d; N_prop = %d, %d\n', nx, ny, N_prop_L, N_prop_R);

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

% Build ny-by-ny unitary matrix phi where the a-th column is the a-th transverse mode
phi = channels.fun_phi(channels.kydx);

% Build input matrix B and output matrix C; the -2i prefactor will be multiplied at the end
% B and C at one pixel outside the scattering region (n=0 and n=nx+1), as dense matrices
B_L = phi(:,channels.L.ind_prop(ind_in_L))*spdiags(reshape(channels.L.sqrt_mu(ind_in_L),[],1), 0, N_in_L, N_in_L);
B_R = phi(:,channels.R.ind_prop(ind_in_R))*spdiags(reshape(channels.R.sqrt_mu(ind_in_R),[],1), 0, N_in_R, N_in_R);
% Note that when the frequency k0dx is complex, sqrt_mu is also complex, and we should take the conjugate transpose of phi but not of sqrt_mu
C_L = spdiags(reshape(channels.L.sqrt_mu(ind_out_L),[],1), 0, N_out_L, N_out_L)*(phi(:,channels.L.ind_prop(ind_out_L))');
C_R = spdiags(reshape(channels.R.sqrt_mu(ind_out_R),[],1), 0, N_out_R, N_out_R)*(phi(:,channels.R.ind_prop(ind_out_R))');

% Retarded Green's function G0 of a semi-infinite homogeneous space, evaluated at the surface (just before the space is terminated)
G0_L = phi*spdiags(exp(1i*channels.L.kxdx(:)), 0, ny, ny)*(phi');
if syst.epsilon_R == syst.epsilon_L
    G0_R = G0_L;
else
    G0_R = phi*spdiags(exp(1i*channels.R.kxdx(:)), 0, ny, ny)*(phi');
end

clear phi

% Finite-difference wave operator for one slice, without the (k0dx)^2*epsilon term
% adds 2 for the laplacian in x direction
A0 = -build_laplacian_1d(ny, syst.yBC) + 2*speye(ny);

k0dx2 = (syst.k0dx)^2;

%% Main computation

% Iterate through by adding slices
% To avoid unnecessary steps, we pick the scheme based on which parts of the S-matrix are requested
if N_in_R==0 && N_out_R==0
    % input and output both on the left; loop from right to left
    G_LL = G0_R;
    for n = nx:-1:1
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_LL = inv(A_nn - G_LL);
    end
    G_LL = (eye(ny) - G0_L * G_LL) \ G0_L;
    smatrix = (-2i)*(C_L * (G_LL * B_L));
elseif N_in_R==0
    % input from the left; loop from right to left
    G_LL = G0_R; G_RL = G0_R;
    for n = nx:-1:1
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_LL = inv(A_nn - G_LL);
        G_RL = G_RL * G_LL;
    end
    G_LL = (eye(ny) - G0_L * G_LL) \ G0_L;
    G_RL = G_RL * G_LL;
    smatrix = (-2i)*[(C_L * G_LL) * B_L; C_R * (G_RL * B_L)]; % also works when N_out_L==0
elseif N_in_L==0 && N_out_L==0
    % input and output both on the right; loop from left to right
    G_RR = G0_L;
    for n = 1:nx
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_RR = inv(A_nn - G_RR);
    end
    G_RR = (eye(ny) - G0_R * G_RR) \ G0_R;
    smatrix = (-2i)*(C_R * (G_RR * B_R));
elseif N_in_L==0
    % input from the right; loop from left to right
    G_RR = G0_L; G_LR = G0_L;
    for n = 1:nx
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_RR = inv(A_nn - G_RR);
        G_LR = G_LR * G_RR;
    end
    G_RR = (eye(ny) - G0_R * G_RR) \ G0_R;
    G_LR = G_LR * G_RR;
    smatrix = (-2i)*[C_L * (G_LR * B_R); (C_R * G_RR) * B_R]; % also works when N_out_R==0
else
    % general case: input from both sides; loop from left to right
    % Initial step: retarded Green's function for an semi-infinite homogeneous space on the left
    G_RR = G0_L; G_RL = G0_L; G_LL = G0_L; G_LR = G0_L;

    % Main loop of RGF: attach the scattering region slice by slice
    for n = 1:nx
        A_nn = A0 - spdiags(k0dx2*(syst.epsilon(:,n)), 0, ny, ny);
        G_RR = inv(A_nn - G_RR);
        G_RL = G_RR * G_RL;
        G_LL = G_LL + G_LR * G_RL;
        G_LR = G_LR * G_RR;
    end

    % Final step: attach homogeneous space on the right
    G_RR = (eye(ny) - G0_R * G_RR) \ G0_R;
    G_RL = G_RR * G_RL;
    G_LL = G_LL + G_LR * G_RL;
    G_LR = G_LR * G_RR;

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
