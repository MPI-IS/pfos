function ode = linearPursuit ()
% LINEARPURSUIT - an ODE derived from a linear pursuit (E5 in Hull et al.)

% linearPursuit.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-25
% Version: 0.1

ode = struct();

ode.name = 'Lin. pursuit';

ode.fun = @(t,x) [x(2); sqrt(1 + x(2)^2)/(25 - t)];
ode.tspan = [0., 20.];
ode.x0 = [0.; 0.];

ode.D = 1;
ode.m = 2;

ode.isStiff = 0;

ode.filename = mfilename('fullpath');
[dirname, fname] = fileparts(ode.filename);
solFilename = [dirname, filesep, 'data', filesep, fname, '.mat'];

if ~exist(solFilename, 'file')
warning('hullOdeProb:solutionMissing', ...
        'continuous solution for this ODE not yet computed');
ode.sol = @(t) NaN(numel(ode.x0), numel(t));
else
solStruct = load(solFilename);

ode.sol = @(t) deval(solStruct.sol, t).';
end

ode.hasAnalyticSol = false;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

end % function