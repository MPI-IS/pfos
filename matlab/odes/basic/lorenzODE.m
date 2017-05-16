function ode = lorenzODE (sigma, rho, beta, x0, tspan)
% LORENZODE - Lorenz' chaotic equations for simplified atmospheric
%             convection
%
% Also see the Wikipedia entry for it:
% http://en.wikipedia.org/wiki/Lorenz_system
%
% Input:
%   sigma - system parameter
%   rho   - system parameter
%   beta  - system parameter
%   x0    - initial values
%   tspan - integration domain

% lorenzODE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-26
% Version: 0.1

ode = struct();

ode.name = 'Lorenz'' system';

ode.fun = @(t,x) [-sigma*x(1) + sigma*x(2); 
                  -x(1)*x(3) + rho*x(1) - x(2);
                  x(1)*x(1) - beta*x(3)];
ode.tspan = tspan;
ode.x0 = x0;

ode.D = 3;
ode.m = 1;

ode.isStiff = 0; 

ode.sol = @(t) NaN(size(t));
ode.hasAnalyticSol = false;
ode.solRefAtEnd = NaN(size(ode.x0));

ode.filename = mfilename('fullpath');

end % function