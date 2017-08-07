function bt = Kutta3 ()
% Kutta3 - Returns the Butcher tableau of Kutta's third-order method
%
% Returns:
%   bt - a ButcherTableau object for Kutta's third-order method

% Kutta3.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0., 0., 0.;
       0.5, 0., 0.;
       -1., 2., 0.];
  b = [(1./6.), (2./3.), (1./6.)];
  c = [0. ; 0.5; 1.];
  
  bt = ButcherTableau(W, b, c);

end % function