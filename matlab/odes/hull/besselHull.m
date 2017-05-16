function ode = besselHull ()
% BESSELHULL - Bessel's ODE with selected parameters from Hull et. al

% besselHull.m
% Author: Michael Schober (mschober@tuebingen.mpg.de)
% Date: 2014-07-24
% Version: 0.1

ode = besselODE (0.5, -1);
ode.filename = mfilename('fullpath');

ode.isStiff = 0;

end % function