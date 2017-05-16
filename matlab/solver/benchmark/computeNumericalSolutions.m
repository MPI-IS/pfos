function computeNumericalSolutions (benchmarks, force)
% COMPUTENUMERICALSOLUTIONS - Recomputes the numerical solutions for a
%                             specific benchmark
%
% This function is not programmed defensively, so do not expect it working
% correctly, if you provide errornous input.
%
% Inputs:
%   benchmark - optional: a cell array containing the names of the
%               benchmarks to update as specified by their directory names.
%               Default: updates all benchmarks
%   force     - optional: boolean indicating whether existing solutions
%               should be overwritten. Default: false

% computeNumericalSolutions.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2014-10-10
% Version: 0.1

  if nargin < 2
    force = false;
  end
  
  if nargin < 1
    benchmarks = {'hull' ...
                , 'mazzia' ...
                , 'hw' ...
                };
  end
  
  if ~iscell(benchmarks)
    benchmarks = { benchmarks };
  end
  
  for i=1:numel(benchmarks)
    odes = eval([benchmarks{i}, 'Benchmark']);
    
    disp(['Computing the numerical solution of ', num2str(numel(odes)), ...
         ' problems in the ', benchmarks{i}, ' benchmark.', ...
         ' This might take a while.']);
    
    for j=1:numel(odes)
      
      disp(['Solving ODE ', num2str(j), ' of ', num2str(numel(odes))]);
      solFile = getSolutionFilename(benchmarks{i}, odes{j});
      
      if odes{j}.hasAnalyticSol
        
        disp('... skipping, because it has an analytical solution');
        continue;
        
      elseif ~exist(solFile, 'file') || force
        
        sol = odeGroundTruth(odes{j});
        save(solFile, 'sol');
        
        disp('... done');
        
      else
        
        disp('... skipping, because solution file already exists');
        
      end % if
      
    end % for odes
  end % for benchmarks

end % function

function filename = getSolutionFilename (benchmark, ode)
% GETSOLUTIONFILENAME - Get the correct filename of the solution file
%                       depending on the benchmark. It is slightly
%                       different for Mazzia then for the other benchmarks.

  [path, odeFilename] = fileparts(ode.filename);
  
  if strcmp(benchmark, 'mazzia')
    filename = [path, filesep, 'data', filesep, ode.name, '.mat'];
  else
    filename = [path, filesep, 'data', filesep, odeFilename, '.mat'];
  end

end % function