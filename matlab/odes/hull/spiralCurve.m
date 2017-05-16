function ode = spiralCurve ()
% SPIRALCURVE - A spiral curve ODE problem

% spiralCurve.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-23
% Version: 0.1

ode = struct();

ode.name = 'Spiral curve';

ode.fun   = @(t, x) (x-t)./(x+t);
ode.tspan = [0., 20.];
ode.x0    = 4.;

ode.D = 1;
ode.m = 1;

ode.isStiff = 0;

% TODO Fix this
% ode.sol = @(theta) 4*exp(pi/2)*exp(-theta);
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