function v = vec (A)
% VEC - vectorizes an input matrix
%
% For a MxN matrix, returns a column vector of size MN by stacking columns
% on top of each other.
%
% Inputs:
%   A - MxN matrix to be put in vector form
%
% Returns:
%   v - column vector of size MN containing the elements of A

% vec.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-03
% Version: 0.1
% Purpose: vectorizes an input matrix

  v = A(:);
  
end % function