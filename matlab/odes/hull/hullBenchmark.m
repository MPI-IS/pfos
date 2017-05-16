function odes = hullBenchmark ()
% HULLBENCHMARK - Returns a cell array containing all ODE problems from the
%                 Hull et al. (1972) benchmark set

% hullBenchmark.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-28
% Version: 0.1

odes = { negExp; % the negative exponential % set A1-A5
         riccati; % a special case of the Riccati equation
         oscillatory; % an oscillatory problem
         logisticCurve; % a logistic curve
         spiralCurve; % a spiral curve
         conflictPop; % the growth of two conflicting populations % B1-B5
         linChemReac; % a linear chemical reaction
         nonLinChemReac; % a nonlinear chemical reaction
         torusSurf; % the integral surface of a torus
         rigidBodyMotion; % Euler equations of motion for a rigid body without external forces
         radActDecay1; % a radioactive decay chain % set C1-C5
         radActDecay2; % another radioactive decay chain
         parabolicPDESmall; % derived from a parabolic differential equation
         parabolicPDEBig; % as C3 except with 51 equations
         fiveBodyProblem; % Five body problem: the motion of 5 outer planets about the sun
         orbitEquation1; % orbit equations % set D1-D5
         orbitEquation2;
         orbitEquation3;
         orbitEquation4;
         orbitEquation5;
         besselHull; % derived from Bessel's equation of order 1/2 with the origin shifted one unit to the left % E1-E5
         vanDerPolHull; % derived from Van der Pol's equation
         duffingHull; % derived from Duffing's equation
         fallingBody; % derived from the falling body equation
         linearPursuit % derived from a linear pursuit equation
         };

end % function