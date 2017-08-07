function bt = Heun ()
% HEUN - Returns the Butcher tableau of Heun's method
%
% Returns:
%   bt - a ButcherTableau object for Heun's method

% Heun.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0. , 0.;
       1. , 0.];
  b = [0.5, 0.5];
  c = [0. ; 1. ];
  
  bt = ButcherTableau(W, b, c);

end % function