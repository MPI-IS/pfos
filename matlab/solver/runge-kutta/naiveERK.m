function [tOut, xOut, YAll, TAll] = naiveERK (odefun, tspan, x0, bt, h)
% NAIVEERK - Implements a naive Runge-Kutta implementation w/o step-size
%            control
%
% Inputs:
%   odefun - a function handle of the form odefun(t0, x0) = x0'
%   tspan  - vector indicating t0 and tFinal
%   x0     - initial value of the solution at time t0
%   bt     - a ButcherTableau object indicating the RK method to use.
%            Currently, this defaults to RKMethods.Heun.
%   h      - step length of one RKstep. Currently, this defaults to h=0.25.
%
% Results:
%   tOut   - vector containing all output times
%   xOut   - solution vector at times t0ut. Each row in xOut corresponds to
%            a row in tOut.
%   YAll   - all intermediate values of each step in a DxSxN dimensional
%            matrix, where D is the dimension of the problem, S is the
%            number of stages of the RK method and N is the number of
%            output times in tOut
%   TAll   - SxN dimensional output matrix of intermediate evaluation
%            times

% naiveRK.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-14
% Version: 0.1

% TODO Support a more similar interface in comparison to ode45

  if (nargin < 4)
    bt = RKMethods.Heun;
  end
  
  if (nargin < 5)
    h = 0.25;
  end
  
  D = numel(x0);
  N = ceil((tspan(end) - tspan(1)) / h);
  
  tOut = NaN(N+1,1);
  xOut = NaN(N+1,D);
  
  YAll = NaN(D, stages(bt), N);
  TAll = NaN(stages(bt), N);
  
  % Start values
  tOut(1)   = tspan(1);
  xOut(1,:) = x0;
  
  % Iteration
  for n = 1:N
    [xOut(n+1,:), YAll(:,:,n), tOut(n+1), TAll(:,n)] = ...
      ERKstep(bt, odefun, tOut(n), xOut(n,:).', h);
  end % for n

end % function