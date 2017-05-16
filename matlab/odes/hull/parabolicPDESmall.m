function ode = parabolicPDESmall ()
% PARABOLICPDESMALL - Derived from a parabolic PDE

% parabolicPDESmall.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = parabolicPDE(10);

ode.name = [ode.name, ' small'];

ode.filename = mfilename('fullpath');

ode.isStiff = 0;

end % function