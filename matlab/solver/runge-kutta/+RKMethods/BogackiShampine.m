function bt = BogackiShampine ()
% BOGACKISHAMPINE - Bogacki-Shampine embedded RK method of order 2/3
%
% Bogacki-Shampine is implemented in Matlab in ode23.
%
% Returns:
%   bt - a ButcherTableau object for the Bogacki-Shampine method

% BogackiShampine.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

% TODO Order seems reversed, actually. Check references again.

  W = [0., 0., 0., 0.;
       0.5, 0., 0., 0.;
       0., 0.75, 0., 0.;
       (2./9.), (1./3.), (4./9.), 0.];
  b = [(2./9.), (1./3.), (4./9.), 0.];
  c = [0.; 0.5; 0.75; 1.];
  
  bStar = [(7./24.), 0.25, (1./3.), (1./8.)];
  
  bt = ButcherTableau(W, b, c, bStar);
  
end % function