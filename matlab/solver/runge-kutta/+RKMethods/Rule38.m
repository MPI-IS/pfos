function bt = Rule38 ()
% RULE38 - Returns the Butcher tableau of the almost-as-famous 3/8th rule
%          by Kutta from his same 1991 paper.
%
% Returns:
%   bt - a ButcherTableau object for the 3/8th rule

% Rule38.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0., 0., 0., 0.;
       (1./3.), 0., 0., 0.;
       -(1./3.), 1., 0., 0.;
       1., -1.,  1., 0.];
  b = [(1./8.), (3./8.), (3./8.), (1./8.)];
  c = [0. ; (1./3.); (2./3.); 1.];
  
  bt = ButcherTableau(W, b, c);

end % function