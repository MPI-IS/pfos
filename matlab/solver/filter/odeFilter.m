function [tout, xout, stats, mf, Pf, model] = ...
  odeFilter (odefun, tspan, x0, options)
% ODEFILTER - an ODE solver based on probabilistic filtering

% odeFilter.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-11-03
% Version: 0.1

%% Configuration

startTime = tic;
odeStartTime = NaN(1); % pre-defined to not allocate in loop

stats = struct ();
stats.fcnCalls = 0;
stats.totalOdeTime = 0;

if nargin < 4
  options = struct ();
end

if isfield(options, 'N')
  N = options.N;
else
  N = 3; % four times IWP is not as stable as three times IWP
end

D = numel(x0);

F     = diag(1:N, 1);
L     = [zeros(N, 1); 1/factorial(N)];
B     = diag(1./factorial(0:N));
Hobs  = [0, 1, zeros(1, N-1)];
Hpred = [1, zeros(1, N)];

Id    = eye(N+1); % used in the update of the covariance matrix

% Error control parameters
if isfield(options, 'ewt')
  ewt = options.ewt;
else
  if isfield(options, 'RelTol')
    rtol = options.RelTol;
  else
    rtol = 1e-3;
  end
  
  if isfield(options, 'AbsTol')
    atol = options.AbsTol;
  else
    atol = 1e-6;
  end
  
  ewt = @(y) 1 ./ (rtol .* abs(y) + atol);
end

% Hyper-parameters
if isfield(options, 'lambda') % learning rate
  lambda = options.lambda;
else
  lambda = 1.;
end

if isfield(options, 'alpha')
  alpha = options.alpha;
else
  alpha = 0.;
end

if isfield(options, 'ErrorPerUnitStep') && options.ErrorPerUnitStep
  ErrorPerUnitStep = true; 
else
  ErrorPerUnitStep = false;
end

if isfield(options, 'eta')
  eta_min = options.eta(1);
  eta_max = options.eta(2);
else
  eta_min = 0.1;
  eta_max = 5.0;
end

if isfield(options, 'rho')
  rho = options.rho;
else
  rho = 0.9;
end

if isfield(options, 'nu')
  nu = options.nu;
else
  nu = 0.75;
end

if isfield(options, 'k_I')
  k_I = options.k_I;
else
  k_I = 0.06 / 2;
end

if isfield(options, 'k_P')
  k_P = options.k_P;
else
  k_P = 0.08 / 2;
end

h_min = 1e-8;

lastcheck = tic();

%% Setup

% Initialize output variables
k = 1;
K = 128;

tout = NaN(K,1); tout(k)   = tspan(1);
xout = NaN(K,D); xout(k,:) = x0.';

% Constants depending on order used
switch N
  case 1
    c = 0;
    W = 0;
  case 2
    u = 2/3;
    c = [0; u];
    W = [0, 0; 
         u, 0];
  case 3
    u = 4/9;
    v = 2/3;
    c = [0; u; v];
    W = [0, 0, 0;
         u, 0, 0;
         (v-(v*(v-u))/(u*(2-3*u))), (v*(v-u))/(u*(2-3*u)), 0];
  case 4
    u = 1/3;
    v = 1/2;
    c = [0.; u; v; 1];
    % In this case, W is a bit more elaborated
    b1 = (1-2*(u+v) + 6*u*v)/(12*u*v);
    b2 = (2*v - 1)/(12*u*(1-u)*(v-u));
    b3 = (1 - 2*u)/(12*v*(1-v)*(v-u));
    b4 = (3-4*(u+v) + 6*u*v)/(12*(1-u)*(1-v));
    
    b = [b1, b2, b3, b4];
    
    W = zeros(4);
    
    W(4,3) = b(3)*(1-c(3))/b(4);
    
    W([3 4],2) = [b(3)*c(2)*c(3), b(4)*c(2)*c(4); b(3), b(4)] \ ...
      [1/8 - W(4,3)*b(4)*c(3)*c(4); b(2)*(1-c(2))];
    
    W(:,1) = c - sum(W(:,2:end),2);
  otherwise
    error('ODEFILTER:illegalArgument', ...
          ['N = ', num2str(N), ' is not supported yet']);
end % switch N

if isfield(options, 'InitialStep')
  [Y, h, m0, P0] = initColDeriv(odefun, tspan(1), x0, c, W, ewt, ...
                                options.InitialStep);
else
  [Y, h, m0, P0] = initColDeriv(odefun, tspan(1), x0, c, W, ewt);
end

% Initialize E_t_old for the PI controller
if ~ErrorPerUnitStep
  E_t_old = rho;
else
  E_t_old = rho * h;
end

% First step has been computed analytically, but in old state-space
m0 = B * m0;
P0 = B * P0 * B.';
P0 = repmat(P0, [1, 1, D]);

if nargout > 3
  mf = NaN(N+1,D,K);
  mf(:,:,k) = m0;
  
  if nargout > 4
    Pf = NaN(N+1,N+1,D,K);
    Pf(:,:,:,k) = P0;
    
    if nargout > 5
      model = struct ();
      
      model.data   = NaN(D,K);
      model.data(:,k) = Y(:,1);
      model.errors = NaN(D,K);
      model.errors(:,k) = Y(:,1);
      model.ssqs   = NaN(D,K);
      model.ssqs(:,k) = ones(D,1);
      
      model.N = N;
      model.alpha  = alpha;
      model.lambda = lambda;
      
      model.ewt = ewt;
      model.ErrorPerUnitStep = ErrorPerUnitStep;
      model.InitialStep = h;
      
      model.h_rec = cell(1,K);
      
    end % nargout > 5 - model defined
    
  end % nargout > 4 - Pf defined
  
end % nargout > 3 - mf defined

m_t = NaN(N+1,D);
P_t = NaN(N+1,N+1);

switch N
  case 1
    m_t = m0;
    P_t = P0(:,:,1); % P0 is already D-dimensional
    
  case 2
    % Analytic mean vector
    m_t(1,:) = (x0 + h * ((1 - 1/(2*u))*Y(:,1) + 1/(2*u)*Y(:,2))).';
    m_t(2,:) = ((1 - 1/u)*Y(:,1) + Y(:,2) ./ u).';
    m_t(3,:) = ((Y(:,2) - Y(:,1))/(h*u)).';
    
    % Analytic covariance matrix
    P_t(1,1) = -h^5*(u^3/24 - u^2/6 + u/6 - 1/20);
    P_t(1,2) = -h^4*(u^3/24 - u^2/4 + u/3 - 1/8); P_t(2,1) = P_t(1,2);
    P_t(1,3) = -h^3*(u^3/24 - u^2/6 + u/3 - 1/6); P_t(3,1) = P_t(1,3);
    P_t(2,2) = (h^3*(u-1)^2)/3;
    P_t(2,3) = h^2*(u^2/6 - (2*u)/3 + 1/2); P_t(3,2) = P_t(2,3);
    P_t(3,3) = -h*((2*u)/3 - 1);
    
  case 3
    % Analytic mean vector
    m_t(1,:) = (x0 + h * ((Y(:,3) - Y(:,1) + 2*v*Y(:,1))/(2*v) ...
              - ((3*v-2)*(Y(:,1) - Y(:,2)))/(6*u*v) ...
              - ((3*v-2)*(Y(:,2) - Y(:,3)))/(6*v*(u-v)))).';

    m_t(2,:) = ((u*(Y(:,1)-Y(:,3)) + v*(Y(:,2)-Y(:,1)) ...
              + u^2*(Y(:,3)-Y(:,1)) + v^2*(Y(:,1)-Y(:,2)) ...
              + Y(:,1)*(u*v*(u-v)))/(u*v*(u-v))).';
            
    m_t(3,:) = ((2*(u*(Y(:,1)-Y(:,3)) + v*(Y(:,2)-Y(:,1))) ...
              + u^2*(Y(:,3)-Y(:,1)) + v^2*(Y(:,1)-Y(:,2)))/(h*u*v*(u-v))).';
            
    m_t(4,:) = (2*(u*(Y(:,1)-Y(:,3)) + v*(Y(:,2)-Y(:,1)))/(h^2*u*v*(u-v))).';
    
    % Analytic covariance matrix
    P_t(1,1) = ...
      (h^7*(- 21*u^4*v^2 + 14*u^4*v - 21*u^3*v^3 + 140*u^3*v^2 - 147*u^3*v ...
      + 42*u^3 - 21*u^2*v^4 + 140*u^2*v^3 - 210*u^2*v^2 + 126*u^2*v ...
      - 28*u^2 - 21*u*v^5 + 140*u*v^4 - 399*u*v^3 + 378*u*v^2 - 112*u*v ...
      + 14*v^5 - 84*v^4 + 210*v^3 - 196*v^2 + 60*v))/(15120*v);
    
    P_t(1,2) = ...
      (h^6*(1 - v)*(u^4*v + u^3*v^2 - 9*u^3*v + 5*u^3 + u^2*v^3 - 9*u^2*v^2 ...
      + 11*u^2*v - 4*u^2 + u*v^4 - 9*u*v^3 + 29*u*v^2 - 16*u*v ...
      - v^4 + 7*v^3 - 18*v^2 + 10*v))/(720*v); P_t(2,1) = P_t(1,2);
    
    P_t(1,3) = ...
      (h^5*(- u^4*v^2 + 2*u^4*v - u^3*v^3 + 8*u^3*v^2 - 18*u^3*v + 8*u^3 ...
      - u^2*v^4 + 8*u^2*v^3 - 24*u^2*v^2 + 24*u^2*v - 8*u^2 - u*v^5 ...
      + 8*u*v^4 - 42*u*v^3 + 72*u*v^2 - 32*u*v + 2*v^5 - 12*v^4 ...
      + 40*v^3 - 56*v^2 + 24*v))/(720*v); P_t(3,1) = P_t(1,3);
    
    P_t(1,4) = ...
      (h^4*(u^4*v + u^3*v^2 - 6*u^3*v + 3*u^3 + u^2*v^3 - 6*u^2*v^2 ...
      + 9*u^2*v - 4*u^2 + u*v^4 - 6*u*v^3 + 27*u*v^2 - 16*u*v + v^5 ...
      - 6*v^4 + 15*v^3 - 28*v^2 + 15*v))/(360*v); P_t(4,1) = P_t(1,4);
    
    P_t(2,2) = (h^5*(u - 1)*(v - 1)^2*(u^2 + u*v + v^2 - 3*v))/(60*v);
    
    P_t(2,3) = ...
      (h^4*(v - 1)*(u^3*v - 3*u^3 + u^2*v^2 - 5*u^2*v + 4*u^2 + u*v^3 ...
      - 11*u*v^2 + 16*u*v - 2*v^3 + 13*v^2 - 15*v))/(120*v);
    P_t(3,2) = P_t(2,3);
    
    P_t(2,4) =  ...
      (h^3*(1 - v)*(u^3 + u^2*v - 2*u^2 + u*v^2 - 8*u*v + v^3 ...
      - 4*v^2 + 10*v))/(60*v); P_t(4,2) = P_t(2,4);
    
    P_t(3,3) = ...
      (h^3*(- u^3*v + 2*u^3 - 2*u^2*v^2 + 6*u^2*v - 4*u^2 - 5*u*v^3 ...
      + 18*u*v^2 - 16*u*v + 10*v^3 - 28*v^2 + 20*v))/(60*v);
    
    P_t(3,4) = ...
      (h^2*(u^3 + 3*u^2*v - 4*u^2 + 9*u*v^2 - 16*u*v + 5*v^3 ...
      - 28*v^2 + 30*v))/(60*v); P_t(4,3) = P_t(3,4);
    
    P_t(4,4) = (h*(- u^2 - 4*u*v - 7*v^2 + 15*v))/(15*v);
    
  case 4
    % Analytic mean vector
    m_t(1,:) = (x0 + h * ( Y(:,1)*(1 - 2*(u+v) + 6*u*v)/(12*u*v) ...
                    + Y(:,2)*(2*v - 1)/(12*u*(v-u)*(1-u)) ...
                    + Y(:,3)*(1 - 2*u)/(12*v*(v-u)*(1-v)) ...
                    + Y(:,4)*(3 - 4*(u+v) + 6*u*v)/(12*(1-u)*(1-v)))).';
  
    
    m_t(2,:) = Y(:,4).';
                  
    m_t(3,:) = ...
      ((Y(:,1)*(1-u)*(v - 1))/(h*u*v) + (Y(:,2)*(1-v))/(h*u*(v-u)*(1-u)) ...
       + (Y(:,3)*(1-u))/(h*v*(v-u)*(v - 1)) + (Y(:,4)*(2-u))/(h*(1-u)) ...
       + Y(:,4)/(h*(1-v))).';
                  
    m_t(4,:) = ((2*Y(:,1)*(u + v - 2))/(h^2*u*v) ...
              + (2*Y(:,2)*(2-v))/(h^2*u*(v-u)*(1-u)) ...
              + (2*Y(:,3)*(2-u))/(h^2*v*(v-u)*(v - 1)) ...
              + (2*Y(:,4)*(3-u-v))/(h^2*(1-u)*(1-v))).';
                  
    m_t(5,:) = (-(6*Y(:,1))/(h^3*u*v) + (6*Y(:,2))/(h^3*u*(v-u)*(1-u)) ...
              + (6*Y(:,3))/(h^3*v*(v-u)*(v - 1)) + (6*Y(:,4))/(h^3*(1-u)*(1-v))).';

    % Analytic covariance matrix
    P_t(1,1) = ...
      (h^9*(6*u^6*v^2 - 3*u^6*v + 6*u^5*v^3 - 27*u^5*v^2 + 20*u^5*v - 4*u^5 ...
      + 6*u^4*v^4 - 27*u^4*v^3 + 28*u^4*v^2 - 12*u^4*v + 2*u^4 ...
      + 6*u^3*v^5 - 27*u^3*v^4 + 28*u^3*v^3 - 12*u^3*v^2 + 2*u^3*v ...
      + 6*u^2*v^6 - 27*u^2*v^5 + 28*u^2*v^4 + 68*u^2*v^3 - 78*u^2*v^2 ...
      + 20*u^2*v - 9*u*v^6 + 38*u*v^5 - 42*u*v^4 - 48*u*v^3 + 70*u*v^2 ...
      - 20*u*v + 3*v^6 - 13*v^5 + 17*v^4 + 5*v^3 - 15*v^2 + 5*v)) ...
      /(725760*v*(1 - u));
    
    P_t(:,2) = 0.; P_t(2,:) = 0.;
    
    P_t(1,3) = ...
      (h^7*(v - 1)*(3*u^6*v + 3*u^5*v^2 - 16*u^5*v + 6*u^5 + 3*u^4*v^3 ...
      - 16*u^4*v^2 + 14*u^4*v - 4*u^4 + 3*u^3*v^4 - 16*u^3*v^3 ...
      + 14*u^3*v^2 - 4*u^3*v + 3*u^2*v^5 - 16*u^2*v^4 ...
      + 14*u^2*v^3 + 76*u^2*v^2 - 40*u^2*v - 6*u*v^5 + 29*u*v^4 ...
      - 24*u*v^3 - 85*u*v^2 + 50*u*v + 3*v^5 - 14*v^4 + 15*v^3 ...
      + 20*v^2 - 15*v))/(120960*v*(u - 1)); P_t(3,1) = P_t(1,3);
    
    P_t(1,4) = ...
      (h^6*(3*u^6*v^2 - 6*u^6*v + 3*u^5*v^3 - 18*u^5*v^2 + 32*u^5*v ...
      - 10*u^5 + 3*u^4*v^4 - 18*u^4*v^3 + 40*u^4*v^2 - 30*u^4*v ...
      + 8*u^4 + 3*u^3*v^5 - 18*u^3*v^4 + 40*u^3*v^3 - 30*u^3*v^2 ...
      + 8*u^3*v + 3*u^2*v^6 - 18*u^2*v^5 + 40*u^2*v^4 + 50*u^2*v^3 ...
      - 192*u^2*v^2 + 80*u^2*v - 9*u*v^6 + 41*u*v^5 - 69*u*v^4 ...
      - 81*u*v^3 + 271*u*v^2 - 117*u*v + 6*v^6 - 28*v^5 + 50*v^4 + 2*v^3 ...
      - 78*v^2 + 39*v))/(60480*v*(u - 1)); P_t(4,1) = P_t(1,4);
    
    P_t(1,5) = ...
      (h^5*(3*u^6*v + 3*u^5*v^2 - 12*u^5*v + 4*u^5 + 3*u^4*v^3 - 12*u^4*v^2 ...
      + 12*u^4*v - 4*u^4 + 3*u^3*v^4 - 12*u^3*v^3 + 12*u^3*v^2 - 4*u^3*v ...
      + 3*u^2*v^5 - 12*u^2*v^4 + 12*u^2*v^3 + 76*u^2*v^2 - 40*u^2*v ...
      + 3*u*v^6 - 12*u*v^5 + 12*u*v^4 + 16*u*v^3 - 140*u*v^2 + 72*u*v ...
      - 3*v^6 + 13*v^5 - 19*v^4 + 5*v^3 + 45*v^2 - 27*v))/(20160*v*(1 - u));
    P_t(5,1) = P_t(1,5);
    
    % Complete 2nd row and 2nd column is zero, which is already set
    
    P_t(3,3) = ...
      (h^5*(v - 1)^2*(u^4 + u^3*v + u^2*v^2 + u*v^3 - 10*u*v - 2*v^3 ...
      + 2*v^2 + 6*v))/(2520*v);
    
    P_t(3,4) = ...
      (h^4*(v - 1)*(u^5*v - 3*u^5 + u^4*v^2 - 5*u^4*v + 4*u^4 + u^3*v^3 ...
      - 5*u^3*v^2 + 4*u^3*v + u^2*v^4 - 5*u^2*v^3 - 16*u^2*v^2 ...
      + 40*u^2*v - 5*u*v^4 + 15*u*v^3 + 37*u*v^2 - 77*u*v + 5*v^4 ...
      - 15*v^3 - 11*v^2 + 33*v))/(2520*v*(u - 1)); P_t(4,3) = P_t(3,4);
    
    P_t(3,5) = ...
      (h^3*(v - 1)*(- u^5 - u^4*v + 2*u^4 - u^3*v^2 + 2*u^3*v - u^2*v^3 ...
      + 2*u^2*v^2 + 20*u^2*v - u*v^4 + 2*u*v^3 + 5*u*v^2 ...
      - 50*u*v + 2*v^4 - 5*v^3 + 25*v))/(840*v*(u - 1));
    P_t(5,3) = P_t(3,5);
    
    P_t(4,4) = ...
      (h^3*(u^5*v - 2*u^5 + 2*u^4*v^2 - 6*u^4*v + 4*u^4 + 2*u^3*v^3 ...
      - 6*u^3*v^2 + 4*u^3*v + 2*u^2*v^4 + 4*u^2*v^3 - 36*u^2*v^2 ...
      + 40*u^2*v + u*v^5 - 12*u*v^4 - 12*u*v^3 + 104*u*v^2 - 96*u*v ...
      - 2*v^5 + 16*v^4 - 8*v^3 - 48*v^2 + 48*v))/(630*v*(1 - u));
    
    P_t(4,5) = ...
      (h^2*(u^5 + 3*u^4*v - 4*u^4 + 3*u^3*v^2 - 4*u^3*v + 3*u^2*v^3 ...
      + 16*u^2*v^2 - 40*u^2*v + 3*u*v^4 + u*v^3 - 95*u*v^2 + 135*u*v + v^5 ...
      - 10*v^4 + 14*v^3 + 54*v^2 - 81*v))/(420*v*(u - 1)); P_t(5,4) = P_t(4,5);
    
    P_t(5,5) = ...
      (h*(u^4 + u^3*v + u^2*v^2 + 10*u^2*v + u*v^3 + 20*u*v^2 - 60*u*v ...
      + v^4 - 5*v^3 - 15*v^2 + 45*v))/(70*v*(1 - u));
    
end % switch N - don't need an otherwise-case, because it's already checked

m_t = B * m_t;
P_t = B * P_t * B.';

P_t = repmat(P_t, [1, 1, D]);

% From NIPS
ssq = abs(m_t(end,:)).';
if nnz(ssq) < numel(ssq)
  if nnz(ssq) > 0
    ssq = max(ssq, min(ssq(ssq > 0)));
  else
    ssq = ones(size(ssq));
  end
end

if true
  ssq = 9 * ones(size(ssq));
end

if N > 1 % for Euler, we don't take an analytic step
  k = k + 1; tout(k) = tout(k-1) + h; xout(k,:) = m_t(1,:);
  
  if nargout > 3
    mf(:,:,k) = m_t;
    
    if nargout > 4
      Pf(:,:,:,k) = P_t;
      
      if nargout > 5
        
        model.data(:,k)   = m_t(2,:);   % these two are imagined values
        model.errors(:,k) = zeros(D,1); % because they don't always apply
        model.ssqs(:,k)   = ssq;
        
      end % model defined
    end % Pf defined
  end % mf defined
  
end % if

%% Main filter loop
while true
  
  h = min(h, tspan(end) - tout(k));
  
  if h <= 0 % reached the end of the integration interval
    break;
  end
  
  % Prepare step loop
  ssq_prev  = ssq;
  numfailed = 0;
  forceStep = false;
  
  % Loop for trying to take one step
  while true
    
    % Prepare step
    Phi11 = sdePhi11Matrix(F, h);
    Phi12 = sdePhi12Matrix(F, L, 1, h);
    
    A_t = Phi11;
    Q_1 = Phi12 * Phi11.';
    
    m_tp = A_t * m_t;
    
    % Evaluate ODE
    odeStartTime = tic;
    dx = odefun(tout(k) + h, m_tp(1,:).');
    stats.totalOdeTime = stats.totalOdeTime + toc(odeStartTime);
    stats.fcnCalls = stats.fcnCalls + 1;
    
    res = dx - m_tp(2,:).';
    
    % Adjust ssq    
    % Maybe turn into a fancy bsxfun/repmat approach, if bottleneck
    exp_errs = NaN(D,1);
    for d=1:D
      exp_errs(d) = Hobs * A_t * P_t(:,:,d) * A_t.' * Hobs.';
    end
    
    ssq_ml = max(res.^2 - alpha * exp_errs, 0.) ./ Q_1(2,2);
    ssq    = lambda * ssq_ml + (1 - lambda) * ssq_prev;
    
    % Error control
    E_t = sqrt(ssq * Q_1(1,1));
    D_t = max(E_t .* ewt(m_tp(1,:).'));
    
    % Either do error control per step or per unit step
    if ~ErrorPerUnitStep
      tol = 1;
    else
      tol = h;
    end
    
    if D_t > tol && ~forceStep % step rejected
      numfailed = numfailed + 1;
      if nargout > 5
        model.h_rec{k} = [model.h_rec{k}, h];
      end
      h = nu * h;
      if h < h_min
        % Replace by more robust mechanism
        warning('Had to resort to extremely small step size h=%.2e',h_min);
        h = min(h_min,tspan(end) - tout(k));
        forceStep = true;
      end
    else
      dtol = rho * tol;
      eta = max(eta_min, min(eta_max, ...
                (dtol/D_t)^(k_I+k_P) * (E_t_old/D_t)^(k_P)));
      E_t_old = D_t;
      break;
    end % if D_t > tol
       
  end % while true - loop to take one step
  
  % Kalman step
  for d=1:D
    P_tp = A_t * P_t(:,:,d) * A_t.' + ssq(d) * Q_1;
    S_tp = Hobs * P_tp * Hobs.';
    K_tp = P_tp * Hobs.' / S_tp;
    
    m_t(:,d) = m_tp(:,d) + K_tp * res(d);
    % update below is numerically more stable than original update
    % P_t(:,:,d) = P_tp - K_tp * S_tp * K_tp.';
    P_t(:,:,d) = (Id - K_tp * Hobs) * P_tp;
    
    % it might be necessary to symmetrize sometimes
    % P_t(:,:,d) = 0.5 * (P_t(:,:,d) + P_t(:,:,d).');
  end % for D
  
  % check whether we need to increase memory
  if k == K
    tout = cat(1, tout, NaN(K, 1));
    xout = cat(1, xout, NaN(K, D));
    
    if nargout > 3
      mf = cat(3, mf, NaN(N+1,D,K));
      
      if nargout > 4
        Pf = cat(4, Pf, NaN(N+1,N+1,D,K));
        
        if nargout > 5
          model.data   = cat(2, model.data, NaN(D,K));
          model.errors = cat(2, model.errors, NaN(D,K));
          model.ssqs   = cat(2, model.ssqs, NaN(D,K));
          model.h_rec  = cat(2, model.h_rec, cell(1,K));
          
        end % model defined
        
      end % Pf defined
      
    end % mf defined
    
    % Unlimited exponential growth is a bad idea, so go linear from a point
    if K < 2^12
      K = K*2;
    else
      K = K+k;
    end
    
  end % if k == K
  
  k = k + 1; tout(k) = tout(k-1) + h; xout(k,:) = m_t(1,:);
  
  if nargout > 3
    mf(:,:,k) = m_t;
    
    if nargout > 4
      Pf(:,:,:,k) = P_t;
      
      if nargout > 5
        model.data(:,k)   = dx;
        model.errors(:,k) = res;
        model.ssqs(:,k)   = ssq;
        
      end % model defined
    end % Pf defined
  end % mf defined
  
  if numfailed == 0.
    h = eta * h;
    if h < h_min
      % Replace by more robust mechanism
      warning('Had to resort to extremely small step size h=%.2e',h_min);
      h = min(h_min,tspan(end) - tout(k));
    end
  end
  
  % progress report
  if toc(lastcheck) > 60
    curperc = (tout(k) - tspan(1)) / (tspan(end) - tspan(1));
    info(sprintf('%2.1f percent solved so far', curperc*100), 1);
    lastcheck = tic();
  end
  
end

%% Finalize

tout = tout(1:k);
xout = xout(1:k,:);

if nargout > 3
  mf = mf(:,:,1:k);
  
  if nargout > 4
    Pf = Pf(:,:,:,1:k);
    
    if nargout > 5
      model.data   = model.data(:,1:k);
      model.errors = model.errors(:,1:k);
      model.ssqs   = model.ssqs(:,1:k);
      model.h_rec  = model.h_rec(:,1:k);
      
    end % model defined
    
  end % Pf defined
  
end % mf defined

stats.totalTime = toc(startTime);
stats.totalOvhd = stats.totalTime - stats.totalOdeTime;

end % function