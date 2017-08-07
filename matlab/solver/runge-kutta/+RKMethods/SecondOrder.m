function bt = SecondOrder(alpha)
% SECONDORDER - Returns the Butcher tableau of any second-order two-stage
%               method as defined by parameter alpha
%
% See Wikipedia for details:
% http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_method#Second-order_methods_with_two_stages
%
% Inputs:
%   alpha - Parameter of the two-stage RK method
%
% Returns:
%   bt - a ButcherTableau object for Ralston's method

% SecondOrder.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  lambda = (1 / (2 * alpha)); % -- ratio parameter in b

  W = [0., 0.;
       alpha, 0.];
  b = [(1 - lambda), lambda];
  c = [0.; alpha];
  
  bt = ButcherTableau(W, b, c);

end % function