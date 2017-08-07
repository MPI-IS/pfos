function bt = ButcherTableau (W, b, c, bStar)
% BUTCHERTABLEAU - Represents a specific Runge-Kutta method via its Butcher
%                  tableau
%
% A Butcher tableau uniquely identifies Runge-Kutta methods. For explicit
% RK methods, it has the following form:
%
% c_1 = 0 | 0
% c_2     | w_21  0
%  *      |  *    *    0
% c_s     | w_s1  *  w_s,s-1 0
% ------------------------------
%  *      | b_1   *   b_s-1  b_s
%
% Currently, the API expects a full s-by-s matrix W and s-dimensional
% vectors c and b. Don't expect the code to check for errors currently.
% 
% Inputs:
%   W  - s-by-s matrix of weights w_{ij}
%   b  - s-dimensional vector of final weights
%   c  - s-dimensional vector of evaluation positions
%
% Returns
%   bt - a Butcher tableau object, containing the weights and evaluation
%        positions of this RK method.

% ButcherTableau.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  bt = struct ();
  
  if (nargin < 4)
    bStar = zeros(size(b));
  end
  
  bt.W = W;
  bt.b = b;
  bt.c = c;
  
  bt.bStar = bStar;
  
  bt.s = numel(bt.b);
  
  bt = class (bt, 'ButcherTableau');

end % function