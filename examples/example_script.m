% A script giving examples of using the cal_smatrix_RGF() function.

%% system parameters of this example
rng default % for reproducibility of random numbers below
ny = 500; nx = 100;  % number of grid points in x and y
syst.k0dx = 2*pi/20;  % k0dx = (2*pi/lambda)*dx where lambda is vacuum wavelength and dx is grid size; use 20 grid points per vacuum wavelength here
syst.epsilon_L = 1.0;  % relative permittivity for the homogeneouse space on the left
syst.epsilon_R = 1.0;  % relative permittivity for the homogeneouse space on the right
syst.epsilon = 1.0 + 2.0*rand(ny, nx);  % relative permittivity of the scattering region; can be anything and can be complex-valued; use random numbers as an example
syst.yBC = 'periodic';  % boundary condition in y

%% calculate the transmission matrix from left to right
in = {'left'};
out = {'right'};
fprintf('computing transmission matrix...\n');
t = cal_smatrix_RGF(syst, out, in);

%% calculate the entire scattering matrix
in = {'left', 'right'};
out = {'left', 'right'};
fprintf('computing full scattering matrix...\n');
[S, channels] = cal_smatrix_RGF(syst, out, in);
N_prop_L = channels.L.N_prop;
N_prop_R = channels.R.N_prop;

% check unitarity of S
err_unitarity = max(max(abs((S'*S) - eye(size(S)))));
fprintf('maximal unitarity violation = %g\n', err_unitarity);

% check symmetry of S, which is guaranteed by reciprocity
% for periodic boundary condition, a permutation is necessary
ind_prop_conj = [channels.L.ind_prop_conj, N_prop_L+channels.R.ind_prop_conj];
S_temp = S(ind_prop_conj,:);
err_reciprocity = max(max(abs(S_temp.' - S_temp)));
fprintf('maximal reciprocity violation = %g\n', err_reciprocity);

% These are the components of S
r = S(1:N_prop_L, 1:N_prop_L); % L->L
t = S(N_prop_L + (1:N_prop_R), 1:N_prop_L); % L->R
tp = S(1:N_prop_L, N_prop_L + (1:N_prop_R)); % R->L
rp = S(N_prop_L + (1:N_prop_R), N_prop_L + (1:N_prop_R)); % R->R

%% calculate selected columns (incident angles) of the transmission and reflection matrices
% find the channel indices for a few incident angles of interest
angles_in_deg = [-45, 0, 45];
kxdx = channels.L.kxdx(channels.L.ind_prop);
kydx = channels.kydx(channels.L.ind_prop);
ind = round(interp1(atan2(kydx, kxdx), 1:N_prop_L, angles_in_deg/180*pi));
% specify those as the input channels from the left
in = []; in.ind_in_L = ind;
out = {'left', 'right'};
fprintf('computing selected columns of the scattering matrix...\n');
S = cal_smatrix_RGF(syst, out, in);
