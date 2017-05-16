function b = isIn(x, lower, upper)
% ISIN - multi-dimensional test whether a variable is between boundaries
%
% IsIn tests whether lower < x < upper pointwise. Singleton dimensions in
% lower and upper are automatically expended, if necessary.
%
% Inputs:
%   x     - points to test
%   lower - lower interval bounds
%   upper - upper interval bounds
%
% Returns:
%   b     - binary mask of points in the interval

% isIn.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-07
% Version: 0.1

  b = bsxfun(@gt, x, lower) & bsxfun(@lt, x, upper);

end % function