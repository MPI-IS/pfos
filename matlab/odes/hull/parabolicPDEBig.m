function ode = parabolicPDEBig ()
% PARABOLICPDEBIG - Derived from a parabolic PDE

% parabolicPDEBig.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = parabolicPDE(51);

ode.name = [ode.name, ' big'];

ode.filename = mfilename('fullpath');

ode.isStiff = 0;

end % function