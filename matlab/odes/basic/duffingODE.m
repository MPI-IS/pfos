function ode = duffingODE (d, a, b, c, w, x0, tspan)
% DUFFINGODE - Duffing's equation for a damped and driven oscillator
%
% Also see the Wikipedia entry for it:
% http://en.wikipedia.org/wiki/Duffing_equation
%
% Input:
%   d - damping delta
%   a - stiffness alpha
%   b - amount of non-linearity in the restoring force
%   c - amplitude of the restoring force gamma
%   w - frequency of the periodic driving force

% duffingODE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-25
% Version: 0.1

ode = struct();

ode.name = 'Duffing''s ODE';

ode.fun = @(t,x) [x(2); c*cos(w*t) - d*x(2) - a*x(1) - b*x(1)^3];
ode.tspan = tspan;
ode.x0 = x0;

ode.D = 1;
ode.m = 2;

ode.isStiff = 0; % from Wikipedia: "Any of the various numeric methods
                 % can be used", i.e., it's not stiff

ode.sol = @(t) NaN(size(t));
ode.hasAnalyticSol = false;
ode.solRefAtEnd = NaN(size(ode.x0));

ode.filename = mfilename('fullpath');

end % function