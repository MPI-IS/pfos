function ode = radActDecay1 ()
% RADACTDECAY1 - A radioactive decay chain ODE problem

% radActDecay1.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

% [-1, 0, 0, ...]
% [ 1,-1, 0, ...]
% [ 0, 1, -1,   ]
% [ ...   -1,  0]
% [ ...    1,  0]
A = diag([-1 * ones(9,1); 0]) + diag(ones(9,1),-1);
x0 = [1; zeros(9,1)];

ode = linearODE(A, x0, [0., 20.]);

ode.name = 'Rad. decay 1';

ode.filename = mfilename('fullpath');

ode.isStiff = 0;

end % function