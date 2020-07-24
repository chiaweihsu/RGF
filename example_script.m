% A script giving examples of using the cal_smatrix_RGF() function.

%% system parameters of this example
rng default % for reproducibility of random numbers
ny = 100; nx = 100;  % system size
syst.kdx = 2*pi/10;  % normalized frequency; 10 grid points per vacuum wavelength
syst.epsilon_L = 1.0;  % relative permittivity on the left
syst.epsilon_R = 1.0;  % relative permittivity on the right
syst.epsilon = 1.0 + 2.0*rand(ny, nx);  % relative permittivity of the scattering region; can be complex-valued
syst.BC = 'periodic';  % boundary condition in y

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

%% calculate selected columns (incident angles) of the transmission and the reflection matrices
% find the channel indices for the desired incident angles
angles_in_deg = [-45, 0, 45];
kxdx = channels.L.kxdx(channels.L.ind_prop);
kydx = channels.kydx(channels.L.ind_prop);
ind = round(interp1(atan2(kydx, kxdx), 1:N_prop_L, angles_in_deg/180*pi));
% specify those as incident channels from the left
in = []; in.ind_in_L = ind;
out = {'left', 'right'};
fprintf('computing selected columns of the scattering matrix...\n');
S = cal_smatrix_RGF(syst, out, in);
