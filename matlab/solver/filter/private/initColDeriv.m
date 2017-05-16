function [Y, h, m0, P0] = initColDeriv (odefun, t0, x0, c, W, ewt, h)
% INITCOLDERIV - Initializes the analytical start by collecting the
%                Runge-Kutta evaluations

% initColDeriv.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-11-12
% Version: 0.1

D = numel(x0);
N = numel(c);

Y = NaN(D, N);

% Starting step size as in Hairer&Wanner, p. 169
Y(:,1) = odefun(t0, x0);
% a)
d0 = max(x0 .* ewt(x0));
d1 = max(Y(:,1) .* ewt(Y(:,1)));
% b)
if min(d0, d1) <= 1e-5
  h0 = 1e-6;
else
  h0 = 0.01 * d0 / d1;
end
% c)
Ye = odefun(t0 + h0, x0 + h0 * Y(:,1));
% d)
d2 = max((Ye - Y(:,1)) .* ewt(Ye - Y(:,1))) / h0;
% e)
if max(d1,d2) <= 1e-15
  h1 = max(1e-6, h0 * 1e-3);
else
  h1 = (0.01 / max(d1,d2))^(1/(N+1));
end
% f)
h1 = min(100 * h0, h1);

if nargin < 7
  h = h1;
else
  if h > h1
    warning('INITCOLDERIV:bigStartingStep', ...
            sprintf('h_init = %.2e > %.2e = h_calc', h, h1));
  end
end % if h defined

c = h * c;
W = h * W;

for j=2:N
  Y(:,j) = odefun(t0 + c(j), x0 + Y(:,1:j-1) * W(j,1:j-1).');
end

% construct m0, P0
m0 = [x0.';
      Y(:,1).';
      zeros(N-1, D)];
    
% TODO This needs to be adjusted, put some big numbers for now
% P0 will be transformed into D-dimensional space outside of this function
P0 = [zeros(2,2), zeros(2, N-1);
      zeros(N-1, 2), 1e2 * ones(N-1)];

end % function