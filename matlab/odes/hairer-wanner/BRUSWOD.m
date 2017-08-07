function ode = BRUSWOD ()
% BRUS - Brusselator reaction without diffusion

% BRUSWOD.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2017-07-12
% Version: 0.1

ode = struct();

ode.name = 'BRUSWOD';

ode.fun   = @(~, x) [1 + x(1,:).^2.*x(2,:) - 4*x(1,:); ...
                     3*x(1,:) - x(1,:).^2*x(2,:)];
ode.tspan = [0., 20];
ode.x0    = [1.5; 3.];

ode.D = 2;
ode.m = ones(1, ode.D);

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
