function dx = numericalDifferentiator(f, x, h)
% NUMERICALDIFFERENTIATOR - computes the numerical derivative of f at x
%
% Inputs:
%   f  - function to differentiate
%   x  - input values
%   h  - offset (optional)
%
% Returns:
%   dx - numerical derivative

% numericalDifferentiator.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-03
% Version: 0.1

  if (nargin < 3)
    h = 1e-5; % least numerical error, found empirically
  end
  
  dx = (f(x+h) - f(x-h)) ./ (2*h);
end
