% A script giving examples of using the cal_smatrix_RGF_3D() function.

%% system parameters of this example
rng default % for reproducibility of random numbers below
nx = 40; ny = 40;  % number of grid points in transverse directions
nz = 20; % number of grid points in axial direction
syst.k0dx = 2*pi/20;  % k0dx = (2*pi/lambda)*dx where lambda is vacuum wavelength and dx is grid size; use 20 grid points per vacuum wavelength here
syst.epsilon_L = 1.0;  % relative permittivity for the homogeneouse space on the left
syst.epsilon_R = 1.0;  % relative permittivity for the homogeneouse space on the right
syst.epsilon = 1.0 + 2.0*rand(nx, ny, nz);  % relative permittivity of the scattering region; can be anything and can be complex-valued; use random numbers as an example
syst.xyBC = 'periodic';  % boundary condition in x and y

%% calculate the transmission matrix from left to right
in = {'left'};
out = {'right'};
fprintf('computing transmission matrix...\n');
t = cal_smatrix_RGF_3D(syst, out, in);

%% calculate the entire scattering matrix
in = {'left', 'right'};
out = {'left', 'right'};
fprintf('computing full scattering matrix...\n');
[S, channels] = cal_smatrix_RGF_3D(syst, out, in);
N_prop_L = channels.L.N_prop;
N_prop_R = channels.R.N_prop;

% check unitarity of S
err_unitarity = max(max(abs((S'*S) - eye(size(S)))));
fprintf('maximal unitarity violation = %g\n', err_unitarity);

% These are the components of S
r = S(1:N_prop_L, 1:N_prop_L); % L->L
t = S(N_prop_L + (1:N_prop_R), 1:N_prop_L); % L->R
tp = S(1:N_prop_L, N_prop_L + (1:N_prop_R)); % R->L
rp = S(N_prop_L + (1:N_prop_R), N_prop_L + (1:N_prop_R)); % R->R
