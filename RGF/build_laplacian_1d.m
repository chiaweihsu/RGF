function A = build_laplacian_1d(n, BC)
%build_laplacian_1d   Returns the finite-difference laplacian operator in 1D.
%   A = build_laplacian_1d(n, BC) returns A as a sparse matrix [(d/dx)^2]*(dx^2)
%   discretized through center difference with uniform grid size dx.
%
%   === Input parameters ===
%   n (non-negative integer scalar);
%       Total number of sites.
%   BC (character vector or string (case insensitive) or scalar number):
%       Boundary condition. For character inputs, the options are:
%         'periodic': E_z(i+n) = E_z(i)
%         'Dirichlet' or 'PEC': E_z(0) = E_z(n+1) = 0
%         'Neumann' or 'PMC': E_z(0) = E_z(1), E_z(n+1) = E_z(n)
%         'DirichletNeumann' or 'PECPMC': E_z(0) = 0, E_z(n+1) = E_z(n)
%         'NeumannDirichlet' or 'PMCPEC': E_z(0) = E_z(1), E_z(n+1) = 0
%       When BC is a scalar number, the Bloch periodic boundary condition is
%       used with f(i+n) = f(i)*exp(1i*kx*a) where kx is the wave
%       number. Here, BC = kx*a = kx*dx*n.
%
%   === Returns ===
%   A (sparse matrix):
%       Matrix representation of [(d/dx)^2]*(dx^2).

% Check input parameters
if ~(isscalar(n) && isnumeric(n) && isreal(n) && (round(n)==n) && n>=0)
    error('The number of sites n must be a non-negative integer scalar');
end
if ~((ischar(BC) && isvector(BC)) || ((isstring(BC) || isnumeric(BC)) && isscalar(BC)))
    error('Boundary condition BC has to be a character vector or string, or numeric scalar (for Bloch periodic BC)');
end

% handle periodic and Bloch periodic boundary conditions
if strcmpi(BC, 'Bloch')
    error('To use Bloch periodic boundary condition, set BC to k*a where k is the Bloch wave number and a is the periodicity');
elseif isnumeric(BC)
    ka = BC;
    BC = 'Bloch';
    if ~isreal(ka)
        warning('k*a = %g + 1i*%g is a complex number.', real(ka), imag(ka));
    end
elseif strcmpi(BC, 'periodic')
    ka = 0;
    BC = 'Bloch';
end

% f = [f(1), ..., f(n)].' 
% first derivative of f
if strcmpi(BC, 'Bloch')
    % f(n+1) = f(1)*exp(1i*ka); f(0) = f(n)*exp(-1i*ka)
    % grad_1*f = df = [df(0.5), ..., df(n-0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1),-exp(-1i*ka)*ones(n,1)], [0,-1,n-1], n, n);
elseif strcmpi(BC, 'Dirichlet') || strcmpi(BC, 'PEC') % PEC on both sides
    % f(0) = f(n+1) = 0
    % grad_1*f = df = [df(0.5), ..., df(n+0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1)], [0,-1], n+1, n);
elseif strcmpi(BC, 'Neumann') || strcmpi(BC, 'PMC') % PMC on both sides
    % df(0.5) = df(n+0.5) = 0
    % grad_1*f = df = [df(1.5), ..., df(n-0.5)].'
    grad_1 = spdiags([ones(n-1,1),-ones(n-1,1)], [1,0], n-1, n);
elseif strcmpi(BC, 'DirichletNeumann') || strcmpi(BC, 'PECPMC') % PEC on the low side, PMC on the high side
    % f(0) = 0; df(n+0.5) = 0
    % grad_1*f = df = [df(0.5), ..., df(n-0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1)], [0,-1], n, n);
elseif strcmpi(BC, 'NeumannDirichlet') || strcmpi(BC, 'PMCPEC') % PMC on the low side, PEC on the high side
    % df(0.5) = 0; f(n+1) = 0
    % grad_1*f = df = [df(1.5), ..., df(n+0.5)].'
    grad_1 = spdiags([ones(n,1),-ones(n,1)], [1,0], n, n);
else
    error('Boundary condition BC = %s is not an eligible option', BC);
end

% grad_2*df = d2f = [d2f(1), ..., d2f(n)].'
grad_2 = -ctranspose(grad_1);

A = grad_2*grad_1;

end
