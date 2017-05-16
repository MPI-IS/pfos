function ode = besselODE (alpha, t0, x0, tspan)
% BESSELODE - Bessel's differential equation

% besselODE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = struct();

ode.name = 'Bessel''s ODE';

if nargin < 4
    tspan = [0.; 20.];
    if nargin < 3
        x0 = [besselj(alpha, -t0); ...
              0.5 * (besselj(alpha-1, -t0) - besselj(alpha+1, -t0))];
    end % if x0
end % if tspan

ode.fun = @(t, x) [x(2); - (x(2)/(t-t0) + (1 - (alpha/(t-t0))^2)*x(1))];
ode.tspan = tspan;
ode.x0 = x0;

ode.D = 1;
ode.m = 2;

ode.isStiff = 0; % actually, I don't know whether Bessel is considered to
                 % be stiff or not

% TODO Check whether this actually holds ...
ode.sol = @(t) [besselj(alpha, vec(t - t0)), ...
                0.5 * (besselj(alpha-1, vec(t-t0)) - besselj(alpha+1, vec(t-t0)))];
ode.hasAnalyticSol = true;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

ode.filename = mfilename('fullpath');

end % function