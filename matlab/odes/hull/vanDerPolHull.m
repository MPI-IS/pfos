function ode = vanDerPolHull ()
% VANDERPOLHULL - Van der Pol's equation with mu=1, x0 = [2; 0] from Hull

% vanDerPolHull.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = vanDerPolODE(1, [2; 0], [0.; 20]);

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