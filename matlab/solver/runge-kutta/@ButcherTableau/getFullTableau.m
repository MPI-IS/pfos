function T = getFullTableau (bt)
% GETFULLTABLEAU - returns the full Butcher tableau (including a possible
%                  embedded b^star) in matrix form
%
% Inputs:
%   bt - a ButcherTableau object
%
% Returns:
%   T  - a matrix [c, W; 0, b.'; 0, b^star.']

% getFullTableau.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  T = [bt.c, bt.W;
       0.,   bt.b];
  
  if (isEmbedded(bt))
    T = [ T;
         0., bt.bStar];
  end % if
  
end % function