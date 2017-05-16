function xsHat = odeSeqGT (ode, ts, xs)
% ODESEQGT- Generates ground truth assuming a fresh IVP for each
%           integration step
%
% Hull evaluates according to local errors, i.e., he collects error
% measurements based on the error of a step, *assuming* that this step was
% a new IVP. Thus, he doesn't penalize errors due to one very bad step.
% However, the error at the end of the integration interval may be
% arbitrarily bad according to this evaluation scheme.
%
% Inputs:
%   ode   - structure containing the original problem description
%   ts    - time knots of the integration steps
%   xs    - estimated solution at times ts
%
% Returns:
%   xsHat - high accuracy solutions for each integration step, assuming an
%           initial value of the previous solution step

% odeSeqGT.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-10-08
% Version: 0.1

  % normalize input
  D = numel(ode.x0);
  if size(xs,2) ~= D
    xs = xs.';
  end

  xsHat = NaN(size(xs));
  xsHat(1,:) = ode.x0.';
  
  numSteps = numel(ts) - 1;
  tspans = [ts(1:numSteps), ts(2:end)];
  
  % parfor i=1:numSteps
  for i=1:numSteps
    odeCopy       = ode;
    odeCopy.tspan = tspans(i,:);
    odeCopy.x0    = xs(i,:).';
    
    sol = odeGroundTruth(odeCopy);
    xsHat(i+1,:) = deval(sol, odeCopy.tspan(end)).';
  end

end % function