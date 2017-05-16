function ode = conflictPop ()
% CONFLICTPOP - ODE model for two conflicting populations

% conflictPop.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-23
% Version: 0.1

ode = struct();

ode.name = 'Confl. Popls.';

ode.fun = @(t, x) [2*(x(1) - x(1)*x(2));
                    -(x(2) - x(1)*x(2))];
ode.tspan = [0., 20.];
ode.x0 = [1.; 3.];

ode.D = 2;
ode.m = [1, 1];

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