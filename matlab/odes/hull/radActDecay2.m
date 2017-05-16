function ode = radActDecay2 ()
% RADACTDECAY2 - Another radioactive decay chain ODE problem

% radActDecay2.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-24
% Version: 0.1

% [-1, 0, 0, ...]
% [ 1,-2, 0, ...]
% [ 0, 2, -3,   ]
% [ ...   -9,  0]
% [ ...    9,  0]
A = diag([-1 * (1:9), 0]) + diag((1:9),-1);
x0 = [1; zeros(9,1)];

ode = linearODE(A, x0, [0., 20.]);

ode.name = 'Rad. decay 2';

ode.filename = mfilename('fullpath');

ode.isStiff = 0;

end % function