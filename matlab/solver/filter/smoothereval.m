function [mt, Pt] = smoothereval (ts, tout, m, P, model)
% SMOOTHEREVAL - Evaluates a posterior distribution from a smoother at
%                points

% smoothereval.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-08-25
% Version: 0.1

error('Code is bugged. Use odeSmoother until this has been fixed');

if ts(1) < tout(1)
  error('Cannot predict in the past for this model');
end

[N, D, ~] = size(m);

F = diag(1:(N-1),1);
L = [zeros(N-1,1); 1/factorial(N-1)];

mt = NaN([N, D, numel(ts)]);
Pt = NaN([N, N, D, numel(ts)]);

K = numel(tout);
kout = 1;

for ks = 1:numel(ts)
  
  kout = find(tout(kout:K) <= ts(ks), 1, 'last') + kout - 1;
  
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
    
    if kout < K
      h = tout(kout+1) - ts(ks);
      
      A = sdePhi11Matrix(F, h);
      Q = sdePhi12Matrix(F, L, 1, h) * A.';
      
      for d=1:D
        Ptp = A * Pt(:,:,d,ks) * A.' + model.ssqs(d,kout) * Q;
        G = Pt(:,:,d,ks) * A.' / Ptp;
        mt(:,d,ks)   = mt(:,d,ks) + G * (m(:,d,kout+1) - A*mt(:,d,ks));
        Pt(:,:,d,ks) = Pt(:,:,d,ks) + G * (P(:,:,d,kout+1) - Ptp) * G.';
      end
    end
    
  end
  
end

end