function [mt, Pt] = filtereval (ts, tout, m, P, model)
% FILTEREVAL - Evaluates a posterior distribution from a filter at points

% filtereval.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-04-30
% Version: 0.1

if ts(1) < tout(1)
  error('Cannot predict in the past for this model');
end

[N, D, ~] = size(m);

F = diag(1:(N-1),1);
L = [zeros(N-1,1); 1/factorial(N-1)];

mt = NaN([N, D, numel(ts)]);
Pt = NaN([N, N, D, numel(ts)]);

kout = 1;

for ks = 1:numel(ts)
  
  kout = find(tout(kout:end) <= ts(ks), 1, 'last') + kout - 1;
  
  if ts(ks) == tout(kout)
    
    mt(:,:,ks)   = m(:,:,kout);
    Pt(:,:,:,ks) = P(:,:,:,kout);
    
  else
    
    h = ts(ks) - tout(kout);
    
    A = sdePhi11Matrix(F, h);
    Q = sdePhi12Matrix(F, L, 1, h) * A.';
    
    for d=1:D
      mt(:,d,ks)   = A * m(:,d,kout);
      Pt(:,:,d,ks) = A * P(:,:,d,kout) * A.' + model.ssqs(d,kout) * Q;
    end
    
  end
  
end

end