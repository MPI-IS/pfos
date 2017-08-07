function odes = hwBenchmark ()
% HWBENCHMARK - Returns a cell array containing all ODE problems from
%               Hairer & Wanner "Solving Ordinary Differential Equations,
%               Vol. I"

% hwBenchmark.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-26
% Version: 0.1

odes = { EULR; % Euler's equations of motion for a rigid body
         AREN; % Arenstorf orbit
         LRNZ; % Saltzman-Lorenz equation
         PLEI; % the Pleiades, a celestial mechanics problem
         ROPE; % the movement of a hanging rope
         BRUS; % Brusselator with diffusion
         BRUSWOD; % Brusselator w/o diffusion
       };

end % function
