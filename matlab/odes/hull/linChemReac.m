function ode = linChemReac ()
% LINCHEMREAC - A linear chemical reaction model

% linChemReac.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-23
% Version: 0.1

A = [-1.,  1.,  0.;
      1., -2.,  1.;
      0.,  1., -1.];

x0 = [2.; 0.; 1.];

ode = linearODE(A, x0, [0., 20.]);

ode.name = 'Lin. chem. reac.';

ode.filename = mfilename('fullpath');

ode.isStiff = 0;

end % function