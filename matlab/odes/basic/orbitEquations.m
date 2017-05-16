function ode = orbitEquations (e)
% ORBITEQUATIONS - 2 dimensional second-order ODE problem
%
% Input:
%   e - the eccentricity of the orbit

% orbitEquations.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = struct();

ode.name = ['Orbit equations, e=', num2str(e)];

ode.fun   = @(t, x) ...
    [x(2); -x(1)./((x(1)^2 + x(3)^2).^(3/2)); ...
     x(4); -x(3)./((x(1)^2 + x(3)^2).^(3/2))];
ode.tspan = [0., 20.];
ode.x0    = [1 - e;
                 0;
                 0;
 sqrt((1+e)/(1-e))];

ode.D = 2;
ode.m = [2, 2];

ode.isStiff = 0;

% TODO There is a solution, but is implicitly given: u - e*sin(u) - t == 0)
ode.sol = @(t) NaN(size(t));
ode.hasAnalyticSol = false;
ode.solRefAtEnd = NaN(size(ode.x0));

ode.filename = mfilename('fullpath');

end % function