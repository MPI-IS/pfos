function ode = PLEI ()
% PLEI - the Pleiades, a celestial mechanics problem

% PLEI.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-26
% Version: 0.1

ode = struct();

ode.name = 'PLEI';

ode.fun = @pleiadesODE;
ode.tspan = [0., 3.];
ode.x0 = initialValues();

ode.D = 14;
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

function dx = pleiadesODE (~, x)

% _Id_entity matrix for easy selection of variables
Id = eye(28);
m  = (1:7)';

dx = NaN(size(x));
dx(1:2:28) = Id(2:2:28,:) * x;

x_dist = bsxfun(@minus, Id(1:4:28,:) * x, (Id(1:4:28,:) * x)');
y_dist = bsxfun(@minus, Id(3:4:28,:) * x, (Id(3:4:28,:) * x)');

% use diag to avoid dividing by zero
R = (x_dist.^2 + y_dist.^2).^(3/2) + diag(ones(7,1));

dx(2:4:28) = sum(bsxfun(@times, m, x_dist) ./ R);
dx(4:4:28) = sum(bsxfun(@times, m, y_dist) ./ R);

end % function

function x0 = initialValues ()

x0 = [ 3.;   % x 1 - first star
       0.;   % x'1
       3.;   % y 1
       0.;   % y'1
       3.;   % x 2 - second star
       0.;   % x'2
      -3.;   % y 2
       0.;   % y'2
      -1.;   % x 3 - third star
       0.;   % x'3
       2.;   % y 3
       0.;   % y'3
      -3.;   % x 4 - fourth star
       0.;   % x'4
       0.;   % y 4
      -1.25; % y'4
       2.;   % x 5 - fifth star
       0.;   % x'5
       0.;   % y 5
       1.;   % y'5
      -2.;   % x 6 - sixth star
       1.75; % x'6
      -4.;   % y 6
       0.;   % y'6
       2.;   % x 7 - seventh star
      -1.5;  % x'7
       4.;   % y 7
       0.];  % y'7
       

end % function