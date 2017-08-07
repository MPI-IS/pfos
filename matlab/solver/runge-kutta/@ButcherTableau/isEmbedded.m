function b = isEmbedded (bt)
% ISEMBEDDED - tests whether there is a second RK embedded in a BT object
%
% Inputs:
%   bt - a ButcherTableau object
%
% Returns:
%   b  - logical 1 (true), iff this BT object contains a vector b^star

% isEmbedded.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-13
% Version: 0.1

  b = any(bt.bStar ~= zeros(size(bt.b)));

end % function