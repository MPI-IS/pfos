function stats = evalSolHull (ode, ts, xs, tau, EPUS)
% EVALSOLHULL - Evaluates the solution to an ODE with the criteria definded
%               by Hull
%
% A deceived step is a step where
%   || xs_i - x(t_i, t_{i-1}, xs_{i-1} ||_inf > tau * |t_i - t_{i-1}|
% and maximum error is
%   || xs_i - x(t_i, t_{i-1}, xs_{i-1} ||_inf / tau
%
% Also see odeSeqGT.m as Hull uses a step-wise error evaluation scheme.
%
% Inputs:
%   ode   - the ode structure of the solved problem
%   ts    - a vector containing the time knots of the solution steps
%   xs    - a matrix containing the solution at the time knots
%   tau   - error tolerance
%   EPUS  - error per unit step (optional): if 1, error units are step
%           sizes, otherwise the error units are just steps
%
% Returns:
%   stats - a cell array containing {1} the local true solution, {2} the
%           number of steps taken, {3} the percent of steps for which the
%           local error exceeded the tolerance, {4} the maximum local
%           truncation error per (unit) step in units of the tolerance and
%           {5} the local truncation errors per (unit) step in units of the
%           tolerance

% evalSolHull.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-10-08
% Version: 0.1

  if nargin < 5
    EPUS = 0.;
  end
  
  stats = cell(5, 1);
  stats{2} = numel(ts) - 1;
  
  xsHat = odeSeqGT(ode, ts, xs);
  stats{1} = xsHat;
  
  if EPUS
    errorUnits = diff(ts);
  else
    errorUnits = ones(numel(ts) - 1, 1);
  end
  
  errors = xs(2:end,:) - xsHat(2:end,:);
  
  errors = max(abs(errors),[],2);
  % errors = sqrt(sum(errors.^2, 2));
  
  errors = bsxfun(@rdivide, errors, errorUnits);
  
  stats{3} = 100 * sum(errors > tau) / stats{2};
  stats{4} = max(errors ./ tau);
  stats{5} = [0; errors];

end % function
