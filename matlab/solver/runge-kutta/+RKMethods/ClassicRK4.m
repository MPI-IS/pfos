function bt = ClassicRK4 ()
% ClassicRK4 - Returns the Butcher tableau of the classical Runge-Kutta
%              method of 4th-order
%
% Returns:
%   bt - a ButcherTableau object for the classical Runge-Kutta method

% ClassicRK4.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0., 0., 0., 0.;
       0.5, 0., 0., 0.;
       0., 0.5, 0., 0.;
       0., 0.,  1., 0.];
  b = [(1./6.), (1./3.), (1./3.), (1./6.)];
  c = [0. ; 0.5; 0.5; 1.];
  
  bt = ButcherTableau(W, b, c);

end % function