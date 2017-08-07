function ode = ROPE ()
% ROPE - the movement of a hanging rope of length 1 under various forces

% ROPE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-27
% Version: 0.1

ode = struct();

ode.name = 'ROPE';

ode.D = 40;
ode.m = 2 * ones(1, ode.D);

ode.fun = @ropeODE;
ode.tspan = [0., 3.723];
ode.x0 = zeros(sum(ode.m),1);

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

function dx = ropeODE (t, x)

F_y = @(t) (1./cosh(4*t - 2.5)).^4;
F_x = 0.4;

n = size(x,1) / 2;
l = (1:n)'; % [1, 2, ..., n]

% x has layout [t_1, t'_1, t_2, t'_2, ..., t_n, t'_n]
theta = x(1:2:end); % [t_1, t_2, t_3, ..., t_n]
dtheta = x(2:2:end); % [t'_1, t'_2, t'_3, ..., t'_n]

% Prepare matrix D and C
d = diff(theta); % [t_2 - t_1, t_3 - t_2, ..., t_n - t_{n-1}]
D = diag(-sin(d), -1) + diag(-sin(-d), 1);

C_d = 2 * ones(n, 1); C_d(1) = 1; C_d(end) = 3;
C = diag(-cos(d), -1) + diag(-cos(-d), 1) + diag(C_d);

dx = NaN(size(x));
dx(1:2:end) = dtheta;

% a)
v = -n * (n+0.5-l) .* sin(theta) - n^2 * F_x .* sin(theta) ...
  + n^2 * cos(theta) * F_y(t) .* double(l <= 3/4 * n);

% b)
w = D * v + dtheta.^2;

% c)
u = C \ w;

% d)
dx(2:2:end) = C * v + D * u;

end % function