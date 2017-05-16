function ode = fiveBodyProblem ()
% FIVEBODYPROBLEM - the motion of 5 outer planets about the sun

% fiveBodyProblem.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-07-25
% Version: 0.1

ode = struct();

ode.name = '5 body prob.';

ode.fun = @fiveBodyProblemODE;
ode.tspan = [0., 20.];
ode.x0 = initialValues();

ode.D = 15;
ode.m = 2 * ones(1, ode.D);

ode.isStiff = 0;

ode.filename = mfilename('fullpath');
[dirname, fname] = fileparts(ode.filename);
solFilename = [dirname, filesep, 'data', filesep, fname, '.mat'];

if ~exist(solFilename, 'file')
warning('hullOdeProb:solutionMissing', ...
        'continuous solution for this ODE not yet computed');
ode.sol = @(t) NaN(numel(ode.x0), numel(t));
else
solStruct = load(solFilename);

ode.sol = @(t) deval(solStruct.sol, t).';
end

ode.hasAnalyticSol = false;
ode.solRefAtEnd = ode.sol(ode.tspan(end));

end % function

function dx = fiveBodyProblemODE (~, x)
% FIVEBODYPROBLEMODE - the actual ODE to solve

% Constant values
g  = 2.95912208286;          % -- gravitational constant
m0 = 1.00000597682;          % -- mass of the sun and the 4 inner planets
m  = [.000954786104043,  ... % -- mass of Jupiter
      .0002855583733151, ... % -- mass of Saturn
      .0000437273164546, ... % -- mass of Uranus
      .0000517759138449, ... % -- mass of Neptune
      .00000277777777778];   % -- mass of Pluto
  
% Little helper such that the extraction (see below) is a one-liner
select = @(A, ridx, cidx) A(ridx,cidx);

% Extracts only odd-indexed elements and puts them in matrix form
% This was done, because 'squeeze' consumed considerable time
y = reshape(select(eye(30),1:2:30,1:30) * x, 3, 5);

dy = zeros(2, 3, 5);
dy(1,:,:) = reshape(select(eye(30),2:2:30,1:30) * x, 3, 5);

r3 = sum(y.^2).^(3/2);

% TODO De-loopify later, if needed
d3 = zeros(5, 5);
for k = 1:5
    for j = 1:5
        d3(k,j) = sum((y(:,k) - y(:,j)).^2)^(3/2);
    end
end % for k

for i=1:3
    for j=1:5
        
        % sum_{k != j) m_k * [ (y_ik - y_ij)/d_jk^3 - y_ik/r_k^3 ]
        sndTerm = m([1:j-1,j+1:5]) * ...
                  ((y(i,[1:j-1,j+1:5]) - y(i,j))./ d3(j,[1:j-1,j+1:5]) ...
                   - y(i,[1:j-1,j+1:5])./r3([1:j-1,j+1:5])).';
        
        dy(2,i,j) = g * ( -(m0 + m(j))*y(i,j)/r3(j) + sndTerm );
    end % for j=planets
end % for i=dimensions

dx = vec(dy);

end % function

function x0 = initialValues ()
% INITIALVALUES - extra function for cleaner code

x0 = [  3.42947415189;   % y 11 -- first planet
        -.557160570446;  % y'11 -- Jupiter
        3.35386959711;   % y 21
         .505696783289;  % y'21
        1.35494901715;   % y 31
         .230578543901;  % y'31
        6.64145542550;   % y 12 -- second planet
        -.415570776342;  % y'12 -- Saturn
        5.97156957878;   % y 22
         .365682722812;  % y'22
        2.18231499728;   % y 32
         .169143213293;  % y'32
       11.2630437207;    % y 13 -- third planet
        -.325325669158;  % y'13 -- Uranus
       14.6952576794;    % y 23
         .189706021964;  % y'23
        6.27960525067;   % y 33
         .0877265322780; % y'33
      -30.1552268759;    % y 14 -- fourth planet
        -.0240476254170; % y'14 -- Neptune
        1.65699966404;   % y 24
        -.287659532608;  % y'24
        1.43785752721;   % y 34
        -.117219543175;  % y'34
      -21.1238353380;    % y 15 -- fifth planet
        -.176860753121;  % y'15 -- Pluto
       28.4465098142;    % y 25
        -.216393453025;  % y'25
       15.3882659679;    % y 35
        -.0148647893090];% y'35

end % function