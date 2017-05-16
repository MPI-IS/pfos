function errCode = testODEDescription (ode)
% TESTODEDESCRIPTION - Checks whether an ODE struct is in a consistent
%                      state and fullfils the current interface
%
% This function captures the interface description of ODEs. If an ODE
% struct returns zero in this function, it is expected to work with all the
% other code in this package.
%
% Currently, this function is only a necessary, not a sufficient condition.
%
% Inputs:
%   ode     - a structure to check
%
%
% Returns:
%   errCode - 0 iff the ODE struct is consistent fullfils the interface, 
%             1 iff the ODE is logically inconsistent
%             2 iff the input is not a valid ODE type

% testODEDescription.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-10-13
% Version: 0.1

  if ~isstruct(ode)
    warning('testODEDescription:notAnODEstruct', ...
            'Input is not a structure');
    errCode = 2;
    return;
  end
  
  ALLFIELDNAMES = ...
    { 'name', 'the name describing this ODE' ...
    ; 'fun', 'the function handle to the ODE in form @(t,x)' ...
    ; 'tspan', ['the integration domain of this ODE as a row-vector ', ...
                'in the form [t0, tEnd]'] ...
    ; 'x0', 'the initial value of this ODE IVP problem' ...
    ; 'D', 'the number of independent dimensions of this ODE' ...
    ; 'm', ['a row vector containing the order of the ODE in each ', ...
            'dimension'] ...
    ; 'isStiff', ['true, if this ODE is considered to be a stiff ', ...
                  'problem, false otherwise'] ...
    ; 'sol', ['a function handle @(t) that returns the reference ', ...
              'solution at times t'] ...
    ; 'hasAnalyticSol', ['true, if this ODE has an analytic solution', ...
                         ', false otherwise'] ...
    ; 'solRefAtEnd', ['a column-vector containing the reference ', ...
                      'solution at the domain end'] ... 
    ; 'filename', ['the filename of the Matlab function which ', ...
                   'generated this ODE struct'] ...
    };
  
  missingFieldIdxs = ~isfield(ode, ALLFIELDNAMES(:,1));
  if any(missingFieldIdxs)
    missingFields = ALLFIELDNAMES(missingFieldIdxs, 1);
    for i=1:numel(missingFields)
      warning('testODEDescription:missingFieldNames', ...
              [' the field ''', missingFields{i}, ''' is missing']);
    end
    warning('testODEDescription:missingFields', ...
            'there are missing fields in this ODE structure');
    errCode = 2;
    return;
  end
  
  if ode.D ~= numel(ode.m)
    warning('testODEDescription:inconsistentDimensionality', ...
            'ODE struct does not fullfil ode.D = numel(ode.m)');
    errCode = 1;
    return;
  end
  
  if sum(ode.m) ~= numel(ode.x0) || sum(ode.m) ~= numel(ode.solRefAtEnd)
    warning('testODEDescription:inconsistentDimensionality', ...
            'ODE struct does not fullfil sum(ode.m) == numel(ode.x0)');
    errCode = 1;
    return;
  end
  
  % more tests to come
  
  errCode = 0;

end % function