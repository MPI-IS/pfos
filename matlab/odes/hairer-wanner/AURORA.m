function ode = AURORA ()
% AURORA - movement of electrical particle in a magnetic field, like the
%          aurora borealis above Polarcirkeln

% AURORA.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2017-01-13
% Version: 0.1

ode = struct();

ode.name = 'AURORA';

ode.D = 2;
ode.m = 2 * ones(1, ode.D);

ode.fun = @auroraODE;
ode.tspan = [0,1];

% constants mentioned in Hairer & Wanner, Eq. 10.18d
u   = 5 * pi / 4;
g   = -0.5;
z0  = 0.314687;
R0  = 0.257453;
r0  = sqrt(R0^2 + z0^2);
Q0  = 1 - (2*g/R0 + R0 / r0^3)^2;
zp0 = sqrt(Q0) * sin(u);
Rp0 = sqrt(Q0) * cos(u);

ode.x0 = [R0; Rp0; z0; zp0];

ode.polar2xyz = @polar2xyz;

ode.isStiff = 0; 

ode.filename = mfilename('fullpath');
[dirname, fname] = fileparts(ode.filename);
solFilename = [dirname, filesep, 'data', filesep, fname, '.mat'];

if ~exist(solFilename, 'file')
warning('hwOdeProb:solutionMissing', ...
        'continuous solution for this ODE not yet computed');
ode.sol = @(t) NaN(numel(ode.x0), numel(t));
else
solStruct = load(solFilename);

ode.sol = @(t) deval(solStruct.sol, t).';
end

ode.hasAnalyticSol = false;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

end % function

function dx = auroraODE(~,x)

g  = -0.5; % gamma

R  = x(1); 
Rp = x(2);
z  = x(3);
zp = x(4);

r  = sqrt(R^2 + z^2);
C  = (2*g/R + R/r^3);

dx = [ Rp; % Rp
       C * (2*g/R^2 + 3 * R^2 / r^5 - 1/r^3);
       zp;
       C * (3*R*z/r^5);
     ];
end

function xyz = polar2xyz (tout, xout)

    g     = -0.5;
    u     = 5*pi/4;
    phi0  = cos(u); % == sin(u)
    
    Rout  = xout(:,1);
    zout  = xout(:,3);
    
    rout  = sqrt(Rout.^2 + zout.^2);
    dpout = (2*g./Rout + Rout./rout.^3) ./ Rout;
    
    pout  = phi0 + cumsum(dpout(1:end-1) .* diff(tout));
    
    xyz = [Rout(1:end-1) .* cos(pout), ...
           Rout(1:end-1) .* sin(pout), ...
           zout(1:end-1)];
end
