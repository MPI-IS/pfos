function St = sampleFilter (ts, tout, m, P, model, num_samples)
% SAMPLEFILTER - Draws samples from the output of odeFilter

% sampleFilter.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-08-25
% Version: 0.1

% TODO Implement such that smoothing distribution might be returned as well

if ts(1) < tout(1)
  error('Cannot sample in the past for this model');
end

[N, D, ~] = size(m);

F = diag(1:(N-1),1);
L = [zeros(N-1,1); 1/factorial(N-1)];

rand_draw = @(P) semichol(P).' * randn([N, num_samples]);

St = NaN([N, D, numel(ts), num_samples]);

kout = find(tout <= ts(end), 1, 'last');

MSGID = 'MATLAB:nearlySingularMatrix';
WSTATE = warning('OFF', MSGID);
lastwarn(''); % reset lastwarn to get an accurate last warning message

% First, compute posterior at ts(end)
ms = m(:,:,end);
Ps = P(:,:,:,end);

if kout < numel(tout)
  
  for k_i = numel(tout)-1:-1:kout
    h = tout(k_i+1) - tout(k_i);
    
    A = sdePhi11Matrix(F, h);
    Q = sdePhi12Matrix(F, L, 1, h) * A.';
    
    for d=1:D
      m_tp = A * m(:,d,k_i);
      P_tp = A * P(:,:,d,k_i) * A.' + model.ssqs(d,k_i) * Q;
      G    = P(:,:,d,k_i) * A.' / P_tp;
      
      ms(:,d)   = m(:,d,k_i) + G * (ms(:,d) - m_tp);
      Ps(:,:,d) = P(:,:,d,k_i) + G * (Ps(:,:,d) - P_tp) * G.';
    end
  end
  
end

% Second, draw sample at the end
if tout(kout) < ts(end)
  
  h = ts(end) - tout(kout);
  
  A = sdePhi11Matrix(F, h);
  Q = sdePhi12Matrix(F, L, 1, h) * A.';
  
  for d=1:D
    ms(:,d)   = A * ms(:,d);
    Ps(:,:,d) = A * Ps(:,:,d) * A.' + model.ssqs(d,kout)*Q;
  end
end
  
for d=1:D
  St(:,d,end,:) = bsxfun(@plus, ms(:,d), rand_draw(Ps(:,:,d)) );
end

% Third, trace it back to the start

for ks=numel(ts)-1:-1:1
  
  klast = kout;
  kout = find(tout(1:kout) <= ts(ks), 1, 'last');
  
  % Update smoother, if necessary
  if kout < klast
    
    for k_i = klast-1:-1:kout
      h = tout(k_i+1) - tout(k_i);
      
      A = sdePhi11Matrix(F, h);
      Q = sdePhi12Matrix(F, L, 1, h) * A.';
      
      for d=1:D
        m_tp = A * m(:,d,k_i);
        P_tp = A * P(:,:,d,k_i) * A.' + model.ssqs(d,k_i) * Q;
        G    = P(:,:,d,k_i) * A.' / P_tp;
        
        ms(:,d)   = m(:,d,k_i) + G * (ms(:,d) - m_tp);
        Ps(:,:,d) = P(:,:,d,k_i) + G * (Ps(:,:,d) - P_tp) * G.';
      end
    end
    
  end
  
  % Draw sample
  mp = m(:,:,kout);
  Pp = P(:,:,:,kout);
    
  if tout(kout) < ts(ks)
    
    h = ts(ks) - tout(kout);
    
    A = sdePhi11Matrix(F, h);
    Q = sdePhi12Matrix(F, L, 1, h) * A.';
    
    for d=1:D
      mp(:,d)   = A * mp(:,d);
      Pp(:,:,d) = A * Pp(:,:,d) * A.' + model.ssqs(d,kout) * Q;
    end
    
  end
    
  h = ts(ks+1) - ts(ks);
  
  A = sdePhi11Matrix(F, h);
  Q = sdePhi12Matrix(F, L, 1, h) * A.';
  
  for d=1:D
    m_tp = A * mp(:,d);
    P_tp = A * Pp(:,:,d) * A.' + model.ssqs(d,kout) * Q;
    G    = Pp(:,:,d) * A.' / P_tp;
    
    St(:,d,ks,:) = ...
      bsxfun(@plus, mp(:,d), ...
        G * squeeze(bsxfun(@minus, St(:,d,ks+1,:), m_tp)) + ...
        rand_draw(Pp(:,:,d) - G * P_tp * G.') );
    
  end
  
end

warning(WSTATE);
[~, LID] = lastwarn();
if strcmp(MSGID, LID)
  warning('odeSmoother:nearlySingularMatrix', ...
          'Some covariance matrices are nearly singular');
end

end