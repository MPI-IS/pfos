function display (bt)
% DISPLAY - displays the Butcher Tableau
%
% Inputs:
%   bt - the ButcherTableau object to dispay

% display.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  disp('Butcher tableau object:');
  disp(getFullTableau(bt));
  disp('Number of stages:');
  disp(bt.s);

end % function