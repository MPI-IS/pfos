function ode = logisticGrowth (r, K, x0, tspan)
% LOGISTICGROWTH - A logistic growth ODE with rate 'r' and capacity 'K'

% logisticGrowth.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-04-28
% Version: 0.1

ode = struct();

ode.name = 'Logistic Growth';

ode.fun = @(t, x) r*x*(1 - x/K);
ode.tspan = tspan;
ode.x0 = x0;

ode.D = numel(x0);
ode.m = ones(1, ode.D);

ode.isStiff = 1;

ode.sol = @(t) K*x0*exp(r*(t-tspan(1)))./(K+x0*(exp(r*(t-tspan(1))) - 1));
ode.hasAnalyticSol = true;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

ode.filename = mfilename('fullpath');

end % function