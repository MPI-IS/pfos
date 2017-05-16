function ode = oscillatory ()
% OSCILLATORY - An oscillatory ODE problem

% oscillatory.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-23
% Version: 0.1

ode = struct();

ode.name = 'Oscillatory ODE';

ode.fun   = @(t, x) x .* cos(t);
ode.tspan = [0., 20.];
ode.x0    = 1.;

ode.D = 1;
ode.m = 1;

ode.isStiff = 0;

ode.sol = @(t) exp(sin(t));
ode.hasAnalyticSol = true;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

ode.filename = mfilename('fullpath');

end % function