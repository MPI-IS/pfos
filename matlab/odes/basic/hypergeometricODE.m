function ode = hypergeometricODE (c, a, b, x0, tspan)
% HYPERGEOMETRICODE - an ordinary second-order differential equation
%
% Every second-order ordinary differential equation with at most three 
% regular singular points can be rewritten as hypergeometric ordinary
% differential equation.
%
% Also see
%   http://mathworld.wolfram.com/HypergeometricDifferentialEquation.html

% hypergeometricODE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-11-07
% Version: 0.1

ode = struct();

ode.name = 'Hypergeometric ODE';

ode.fun = @(t, x) [x(2); (a*b*x(1) - (c - (a+b+1)*t)*x(2))/(t*(1-t))];
ode.tspan = tspan;
ode.x0 = x0;

ode.D = 1;
ode.m = 2;

% Everything below here are best guesses, based on my current knowledge
ode.isStiff = 1;
ode.sol = @(t) NaN(numel(t),2);
ode.hasAnalyticSol = false;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

ode.filename = mfilename('fullpath');

end % function