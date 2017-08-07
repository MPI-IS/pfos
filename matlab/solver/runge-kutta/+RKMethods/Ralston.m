function bt = Ralston()
% RALSTON - Returns the Butcher tableau of Ralston's method
%
% Returns:
%   bt - a ButcherTableau object for Ralston's method

% Ralston.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  W = [0. , 0.;
       (2./3.), 0.];
  b = [0.25, 0.75];
  c = [0. ; (2./3.)];
  
  bt = ButcherTableau(W, b, c);

end % function