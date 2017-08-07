function bt = HeunEuler ()
% HEUNEULER - Simplest embedded RK method
%
% Returns:
%   bt - a ButcherTableau object for the embedded Heun-Euler method

% HeunEuler.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0. , 0.;
       1. , 0.];
  b = [0.5, 0.5];
  c = [0. ; 1. ];
  
  bStar = [1., 0.];
  
  bt = ButcherTableau(W, b, c, bStar);
  
end % function