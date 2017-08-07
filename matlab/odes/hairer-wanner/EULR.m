function ode = EULR ()
% EULR - Euler equations of motion for a rigid body as given in H&W
%
% Also compare to rigidBodyMotion.m

% EULR.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-26
% Version: 0.1

ode = struct();

ode.name = 'EULR';

f = @(t) double(3*pi <= t & t <= 4*pi) .* (0.25 * (sin(t)).^2);
I = [0.5; 2; 3];
ode.x0 = [1; 0.; 0.9];
ode.tspan = [0., 20.];

ode.fun = @(t, x) [(I(2)-I(3))/I(1)*x(2)*x(3); 
                   (I(3)-I(1))/I(2)*x(3)*x(1);
                   (I(1)-I(2))/I(3)*x(1)*x(2) + f(t)];

ode.D = 3;
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