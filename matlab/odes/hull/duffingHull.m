function ode = duffingHull ()
% DUFFINGHULL - Duffing's equation in Hull's benchmark set
%
% Also see the Wikipedia entry for it:
% http://en.wikipedia.org/wiki/Duffing_equation
%
% This file doesn't use duffingODE, because in Hull's work, sine instead of
% cosine is used - and I don't want to add stupid calculus mistakes.

% duffingHull.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-25
% Version: 0.1

ode = struct();

ode.name = 'Duffing''s ODE';

ode.fun = @(t,x) [x(2); x(1)^3/6 - x(1) + 2*sin(2.78535*t)];
ode.tspan = [0., 20];
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