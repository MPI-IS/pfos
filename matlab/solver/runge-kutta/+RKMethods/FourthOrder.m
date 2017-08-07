function bt = FourthOrder(u, v)
% FOURTHORDER - Returns the Butcher tableau of the fourth order methods for
%               which u ~= v.
%
% See Hairer & Wanner, pt. I, p. 138 for details
%
% Inputs:
%   u  - corresponds to c_2
%   v  - corresponds to c_3
%
% Returns:
%   bt - a ButcherTableau object for Ralston's method

% FourthOrder.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-11-17
% Version: 0.1

  c = [0.; u; v; 1];
  
  b1 = (1-2*(u+v) + 6*u*v)/(12*u*v);
  b2 = (2*v - 1)/(12*u*(1-u)*(v-u));
  b3 = (1 - 2*u)/(12*v*(1-v)*(v-u));
  b4 = (3-4*(u+v) + 6*u*v)/(12*(1-u)*(1-v));
  
  b = [b1, b2, b3, b4];
  
  W = zeros(4);
  
  W(4,3) = b(3)*(1-c(3))/b(4);
  
  W([3 4],2) = [b(3)*c(2)*c(3), b(4)*c(2)*c(4); b(3), b(4)] \ ...
               [1/8 - W(4,3)*b(4)*c(3)*c(4); b(2)*(1-c(2))];
           
  W(:,1) = c - sum(W(:,2:end),2);
  
  bt = ButcherTableau(W, b, c);

end % function