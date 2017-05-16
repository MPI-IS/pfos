function argOut = getOutput (num_output, func, varargin)
% GETOUTPUT - returns a fixed output argument from a multi-output function
% 
% Inputs:
%   num_output - the number of the output argument
%   func       - the function of which the output is extracted
%   varargin   - arguments to the function call
%
% Results:
%   argOut     - the num_output-th output argument of func(varargin)

% getOutput.m
% Author: Michael Schober (michael.schober@tuebingen.mpg.de)
% Date: 2014-03-04
% Version: 0.1

  varargout = cell(num_output,1);
  [varargout{:}] = func(varargin{:});
  argOut = varargout{num_output};

end