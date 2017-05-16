function info (str, lvl, MAX)
% INFO - Display progress information
%
% Inputs:
%   str - formatted string of information to displayed
%   lvl - importance level. If lvl > MAX, it will not be displayed.
%   MAX - Maximum importance level. Default: 1.

% info.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-11-27
% Version: 0.1

persistent msg;

if nargin < 3
  MAX = 1;
  if nargin < 2
    lvl = 0;
    
    if nargin < 1
      fprintf(msg);
      return;
    end
  end
end

msg = sprintf('%s: %s%s\n', datestr(now), repmat('  ',1,lvl), str);

if lvl <= MAX
  fprintf(msg);
end

end