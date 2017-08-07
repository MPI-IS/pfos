function bt = ThirdOrder(u, v)
% THIRDORDER - Returns the Butcher tableau of any third-order three-stage 
%              method as defined by parameters u and v
%
% See Hairer & Wanner, pt. I, p. 142 for details
%
% Inputs:
%   u  - parameter one for third-order methods
%   v  - parameter two for third-order methods
%
% Returns:
%   bt - a ButcherTableau object for Ralston's method

% ThirdOrder.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-04-03
% Version: 0.1

  a32 = v * (v-u) / (u * (2-3*u));
  a31 = v - a32;
  
  b2  = (2-3*v) / (6 * u * (u-v));
  b3  = (2-3*u) / (6 * v * (v-u));
  b1  = 1 - b2 - b3;

  W = [0.,  0.,  0.;
       u,   0.,  0.;
       a31, a32, 0.];
  b = [b1,  b2,  b3];
  c = [0.;  u;   v];
  
  bt = ButcherTableau(W, b, c);

end % function