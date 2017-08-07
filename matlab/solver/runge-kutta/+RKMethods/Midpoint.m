function bt = Midpoint ()
% MIDPOINT - Returns the Butcher tableau of the midpoint method
%
% Returns:
%   bt - a ButcherTableau object for the midpoint method

% Midpoint.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0. , 0.;
       0.5, 0.];
  b = [0. , 1.];
  c = [0. ; 0.5];
  
  bt = ButcherTableau(W, b, c);

end % function