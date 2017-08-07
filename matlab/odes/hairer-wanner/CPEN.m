function ode = CPEN ()
% CPEN - the nonlinear Coupled PENdulum

% CPEN.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-03-23
% Version: 0.1

ode = struct();

ode.name = 'CPEN';

ode.x0 = zeros(4,1);
ode.tspan = [0., 496.];

ode.fun = @coupledPendulumODE;

ode.D = 2;
ode.m = 2 * ones(1, ode.D);

ode.isStiff = 0; 

ode.filename = mfilename('fullpath');
[dirname, fname] = fileparts(ode.filename);
solFilename = [dirname, filesep, 'data', filesep, fname, '.mat'];

if ~exist(solFilename, 'file')
warning('hwOdeProb:solutionMissing', ...
        'continuous solution for this ODE not yet computed');
ode.sol = @(t) NaN(numel(ode.x0), numel(t));
else
solStruct = load(solFilename);

ode.sol = @(t) deval(solStruct.sol, t).';
end

ode.hasAnalyticSol = false;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

end % function

function dx = coupledPendulumODE(t, x)

l1 = 1.;
l2 = 1.;

m1 = 1.;
m2 = 0.99;

rsq = 0.1^2;
c0 = 0.01;

f = @(t) double(abs(t - 1) <= 1) .* sqrt(1 - (1 - t)^2);

dx = NaN(4,1);
dx([1, 3]) = x([2, 4]);

sin1 = sin(x(1));
sin2 = sin(x(3));
cos1 = cos(x(1));
cos2 = cos(x(3));

dx(2) = - sin1 / l1 - c0 * rsq * (sin1 - sin2) * cos1 / m1 / l1^2 + f(t);
dx(4) = - sin2 / l2 - c0 * rsq * (sin2 - sin1) * cos2 / m2 / l2^2;

end % function