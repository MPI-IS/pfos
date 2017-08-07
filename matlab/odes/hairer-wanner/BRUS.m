function ode = BRUS ()
% BRUS - Brusselator reaction-diffusion equation

% BRUS.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-28
% Version: 0.1

ode = struct();

ode.name = 'BRUS';

ode.fun   = @brusODE;
ode.tspan = [0., 7.5];
ode.x0    = initialValues();

ode.D = size(ode.x0, 1);
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

function dx = brusODE(~, x)

N = sqrt(size(x, 1) / 2);
alpha = 2*1e-3;

U = reshape(x(1:N^2), N, N);
V = reshape(x(N^2+1:end), N, N);

% Eq (10.14)
% First, all element-wise operations that can be done w/o extending BCs
Ut = 1 + U.^2 .* V - 4.4 * U - 4 * alpha * (N-1)^2 * U;
Vt = 3.4 * U - U.^2 .* V - 4 * alpha * (N-1)^2 * V;

% Prepare extra boundary condition variables
Ue = [U(:,2), U, U(:,N-1)];
Ue = [Ue(2,:); Ue; Ue(N-1,:)];

Ve = [V(:,2), V, V(:,N-1)];
Ve = [Ve(2,:); Ve; Ve(N-1,:)];

% Now, we can do the rest of the right bracket
Ut = Ut + alpha * (N-1)^2 * (Ue(3:N+2,2:N+1) + Ue(1:N,2:N+1) ...
                           + Ue(2:N+1,3:N+2) + Ue(2:N+1,1:N));
                         
Vt = Vt + alpha * (N-1)^2 * (Ve(3:N+2,2:N+1) + Ve(1:N,2:N+1) ...
                           + Ve(2:N+1,3:N+2) + Ve(2:N+1,1:N));
                         
dx = [vec(Ut); vec(Vt)];

end % function

function x0 = initialValues ()

N = 21;

i = (1:N)'; % row index
j = 1:N;    % col index

% Eq. (10.12)
U = repmat((j-1)./(N-1) + 0.5, N, 1);
V = repmat(5*(i-1)./(N-1) + 1, 1, N);

x0 = [vec(U); vec(V)];

end % function