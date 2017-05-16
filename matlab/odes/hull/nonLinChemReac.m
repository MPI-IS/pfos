function ode = nonLinChemReac ()
% NONLINCHEMREAC - A nonlinear chemical reaction model

% nonLinChemReac.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-23
% Version: 0.1

ode = struct();

ode.name = 'Nonl. chem. reac.';

ode.fun = @(t, x) [-x(1); 
                    x(1) - x(2)^2;
                    x(2)^2];
ode.tspan = [0., 20.];
ode.x0 = [1.; 0.; 0.];

ode.D = 3;
ode.m = ones(1, ode.D);

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