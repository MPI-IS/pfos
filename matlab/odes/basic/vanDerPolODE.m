function ode = vanDerPolODE (mu, x0, tspan)
% VANDERPOLODE - Van der Pol's equation for a non-conservative oscillator

% vanDerPolODE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = struct();

ode.name = 'Van der Pol''s ODE';

ode.fun = @(t,x) [x(2); mu*(1-x(1)^2)*x(2) - x(1)];
ode.tspan = tspan;
ode.x0 = x0;

ode.D = 1;
ode.m = 2;

ode.isStiff = mu > 1e1;

ode.sol = @(t) NaN(size(t));
ode.hasAnalyticSol = false;
ode.solRefAtEnd = NaN(size(x0));

ode.filename = mfilename('fullpath');

end % function