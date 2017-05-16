function sol = odeGroundTruth (ode)
% ODEGROUNDTRUTH - Generates a high accuracy solution of a given ODE

% odeGroundTruth.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-10-08
% Version: 0.1

  % Hairer & Wanner use rtol = atol = h0 in their test cases
  % tol = 1e-15;
  tol = 1e-13;

  matSolOpts = odeset('RelTol', tol ...
                     ,'AbsTol', tol ...
                     ,'MaxStep', 20); % ... % for problems in the Hull data set
                     % ,'MaxOrder', 2 ... % order 2 methods are A-stable
                     % );
                   
  if ode.isStiff
    sol = ode45(ode.fun, ode.tspan, ode.x0, matSolOpts);
  else             
    sol = ode45(ode.fun, ode.tspan, ode.x0, matSolOpts);
  end

end % function