function ode = linearODE (A, x0, tspan)
% LINEARODE - A constructor for linear problems of independent dimensions

% linearODE.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = struct();

ode.name = 'Linear ODE';

ode.fun = @(t, x) A*x;
ode.tspan = tspan;
ode.x0 = x0;

ode.D = numel(x0);
ode.m = ones(1, ode.D);

if numel(A) == 1
  ode.isStiff = A <= -1;
else
  absEV = abs(eig(A));
  ode.isStiff = max(absEV) / min(absEV) > 5.;
end

ode.sol = @(t) trueSol(t, A, x0);
ode.hasAnalyticSol = true;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

ode.filename = mfilename('fullpath');

end % function

function xt = trueSol(t, A, x0)
% TRUESOL - Returns the true solution of the linear ODE system

  xt = NaN(numel(t), numel(x0));
  
  for i=1:numel(t)
      xt(i,:) = expm(A*t(i)) * x0;
  end
end