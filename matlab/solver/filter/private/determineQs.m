function qs = determineQs (F)
% DETERMINEQS - Determines the different number of states for stacked F

% determineQs.m
% Author: Michael Schober
% Date: 2015-06-02
% Version: 0.1

  diagF = diag(F,1);
  qs = diff([0; find(diagF == 0); numel(diagF) + 1]) - 1;

end