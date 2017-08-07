function ode = AREN ()
% AREN - the Arenstorf orbit for the restricted three body problem

% AREN.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-26
% Version: 0.1

ode = struct();

ode.name = 'AREN';

ode.x0 = [0.994; 0.; 0.; -2.00158510637908252240537862224];
ode.tspan = [0., 17.0652165601579625588917206249];

ode.fun = @restrictedThreeBodyProblemODE;

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

function dx = restrictedThreeBodyProblemODE(~, x)

dx = NaN(4,1);
dx([1, 3]) = x([2, 4]);

mu  = 0.012277471;
mut = 1 - mu;

D = [((x(1) + mu)^2  + x(3)^2)^(3/2);
     ((x(1) - mut)^2 + x(3)^2)^(3/2)];
   
dx(2) = x(1) + 2*x(4) - mut*(x(1) + mu)/D(1) - mu*(x(1) - mut)/D(2);
dx(4) = x(3) - 2*x(2) - mut*x(3)/D(1) - mu*x(3)/D(2);

end % function
