function ode = logisticCurve ()
% LOGISTICCURVE - A logistic curve ODE problem

% logisticCurve.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-23
% Version: 0.1

ode = struct();

ode.name = 'Log. curve';

ode.fun   = @(t, x) x/4*(1 - x/20);
ode.tspan = [0., 20.];
ode.x0    = 1.;

ode.D = 1;
ode.m = 1;

ode.isStiff = 0;

ode.sol = @(t) 20./(1 + 19 * exp(-t./4));
ode.hasAnalyticSol = true;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

ode.filename = mfilename('fullpath');

end % function