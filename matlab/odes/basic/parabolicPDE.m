function ode = parabolicPDE (dim, x0, tspan)
% PARABOLICPDE - Derived from a parabolic PDE

% parabolicPDE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

% [-2,  1,  0, ...]
% [ 1, -2,  1, ...]
% [ 0,  1, -2, ...]
% [ ...    -2,   1]
% [ ...     1,  -2]

A = diag(-2 * ones(dim,1)) ...
    + diag(ones(dim-1,1),-1) + diag(ones(dim-1,1),1);

if nargin < 3
    tspan = [0., 20.];
    if nargin < 2
        x0 = [1; zeros(dim-1, 1)];
    end % if ~x0
end % if ~tspan

ode = linearODE(A, x0, tspan);

ode.name = 'Parabolic PDE';

ode.filename = mfilename('fullpath');

end % function