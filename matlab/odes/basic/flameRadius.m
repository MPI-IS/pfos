function ode = flameRadius (delta)
% FLAMERADIUS - An ODE describing the radius of the flame of a match
%
% Taken from the Matlab newsletter at
% http://www.mathworks.de/company/newsletters/articles/stiff-differential-equations.html

% flameRadius.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-10-13
% Version: 0.1

ode = struct();

ode.name = 'Flame radius';

ode.fun   = @(t, x) x.^2 - x.^3;
ode.tspan = [0, 2/delta];
ode.x0    = delta;

ode.D = 1;
ode.m = 1;

ode.isStiff = delta < 1e-4;

ode.sol = @(t) NaN(size(t));
ode.hasAnalyticSol = false;
ode.solRefAtEnd = 1.;

ode.filename = mfilename('fullpath');

end % function