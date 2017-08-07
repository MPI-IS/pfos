function bt = Euler ()
% EULER - Returns the Butcher tableau of the explicit Euler method
%
% Returns:
%   bt - a ButcherTableau object for the explicit Euler method

% Euler.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  bt = ButcherTableau(0., 1., 0.);
  
end % function