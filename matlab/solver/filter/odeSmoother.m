function [ms, Ps, samples] = odeSmoother (tout, m, P, model, num_samples)
% ODESMOOTHER - RTS smoother for output from odeFilter

% odeSmoother.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-04-30
% Version: 0.1

t = numel(tout) - 1;

[N, D, ~] = size(m);

F = diag(1:N-1,1);
L = [zeros(N-1,1); 1/factorial(N-1)];

ms = NaN(size(m));
Ps = NaN(size(P));

ms(:,:,t+1)   = m(:,:,t+1);
Ps(:,:,:,t+1) = P(:,:,:,t+1);

draw_samples = false;

if nargin > 4
  draw_samples = true;
  rand_draw = @(P) semichol(P).' * randn([N, num_samples]);
  
  samples = NaN([size(m),num_samples]);
  for d=1:D
    samples(:,d,t+1,:) = ...
      bsxfun(@plus, m(:,d,t+1), rand_draw(P(:,:,d,t+1)) );
  end
end % if nargin > 4

% TODO Only show RCOND warning once per smoothing
MSGID = 'MATLAB:nearlySingularMatrix';
WSTATE = warning('OFF', MSGID);
lastwarn(''); % reset lastwarn to get an accurate last warning message

while t > 0
  h = tout(t+1) - tout(t);
  
  A_t = sdePhi11Matrix(F, h);
  Q_t = sdePhi12Matrix(F, L, 1, h) * A_t.';
  
  for d=1:D
    m_tp = A_t * m(:,d,t);
    P_tp = A_t * P(:,:,d,t) * A_t.' + model.ssqs(d,t) * Q_t;
    G_t  = P(:,:,d,t) * A_t.' / P_tp;
    
    if draw_samples
      samples(:,d,t,:) = ...
        bsxfun(@plus, m(:,d,t), ...
        G_t * squeeze(bsxfun(@minus, samples(:,d,t+1,:), m_tp)) + ...
        rand_draw(P(:,:,d,t) - G_t * P_tp * G_t.') );
    end
    
    ms(:,d,t)   = m(:,d,t) + G_t * (ms(:,d,t+1) - m_tp);
    Ps(:,:,d,t) = P(:,:,d,t) + G_t * (Ps(:,:,d,t+1) - P_tp) * G_t.';
  end
  
  t = t-1;
end

warning(WSTATE);
[~, LID] = lastwarn();
if strcmp(MSGID, LID)
  warning('odeSmoother:nearlySingularMatrix', ...
          'Some covariance matrices are nearly singular');
end

end