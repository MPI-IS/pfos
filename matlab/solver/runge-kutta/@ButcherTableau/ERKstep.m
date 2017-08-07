function [x1, Y, t1, T, x1Star, err] = ERKstep (bt, odefun, t0, x0, h)
% RKSTEP - Calculates one Runge-Kutta step of step length h with initial
%          values (t0, x0) of the method described by the Butcher tableau
%
% Inputs:
%   bt     - a ButcherTableau object specifying the RK method
%   odefun - a function handle of the form odefun(t0, x0) = x0'
%   t0     - current time
%   x0     - current position
%   h      - step length
%
% Returns:
%   x1     - final prediction of the RK method
%   Y      - s-dimensional vector of intermediate values
%   t1     - next time t1 = t0 + h
%   T      - intermediate time steps of Y
%   x1Star - the embedded solution, or NaN, iff isEmbedded(bt) == 0
%   err    - estimated error, or NaN, iff isEmbedded(bt) == 0

% RKstep.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

% TODO This still needs extensive testing for multi-d!!!

  D = numel(x0);
  T = t0 + h * bt.c;
  Y = NaN(D, bt.s);
  
  for j = 1:bt.s
    Y(:,j) = odefun(T(j), x0 + h * Y(:,1:j-1) * bt.W(j,1:j-1).');
  end % for j
  
  x1 = x0 + h * Y * bt.b.';
  t1 = t0 + h;
  
  if (nargout > 4)
    if (isEmbedded(bt))
      x1Star = x0 + h * bt.bStar * Y;
      err    = abs(x1 - x1Star);
    else
      x1Star = NaN;
      err    = NaN;
    end % if isEmbedded ...
  end % if nargout ...

end % function