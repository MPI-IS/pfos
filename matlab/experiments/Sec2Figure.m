% Sec2Figure.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2015-10-28
% Version: 0.1
% Purpose: Produces the figure for Section 2 in the SIAM paper

%% Setup

clear;
close all;

mpg = [0,0.4717,0.4604]; % color [0,125,122]
dre = [0.4906,0,0]; % color [130,0,0]
ora = [255,153,51] ./ 255;
blu = [0,0,0.509];
gra = 0.5 * ones(3,1);

lightmpg = [1,1,1] - 0.5 * ([1,1,1] - mpg);
lightdre = [1,1,1] - 0.5 * ([1,1,1] - dre);
lightblu = [1,1,1] - 0.5 * ([1,1,1] - blu);
lightora = [1,1,1] - 0.5 * ([1,1,1] - ora);

load handel;
beep = @() sound(y(1:2e4),Fs);
clear y Fs; % -- this is optional, if you don't want a cluttered workspace

%% Configuration

rng(4);

% Use a logistic growth function!
% see
% http://www.scholarpedia.org/article/Linear_multistep_method
% for obvious inspiration

% ODE config

r = 3;
K = 1.;

y0 = 0.1;

t0 = 0.;
T  = 1.5;

tspan = [t0, T];
h = 0.3; N = ceil((T - t0)/h);

ode = @(t,y) r*y.*(1 - y./K);

y   = @(t) K*y0*exp(r*t) ./ (K + y0*(exp(r*t) - 1));
dy  = @(t) K*y0*r*exp(r*t)*(K-y0) ./ (K + y0*(exp(r*t) - 1)).^2; % checked
ddy = @(t) numericalDifferentiator(dy, t); % lazy ;-)

funcs = {y, dy, ddy};

% Solver config

opts = struct ();
opts.lambda = 0.; % don't adjust ssq
opts.ewt = @(y) zeros(size(y)); % always have weighted error < tol
opts.eta = [1, 1]; % scaling of h from one step to another is constant
opts.InitialStep = h;

opts.N = 2;

% Plotting config

num_samples = 3;

b = 0.25;
t = linspace(t0 -b, T +b, 2*261);

fig_y = figure;
fig_dy = figure;
fig_d2y = figure;

figs = {fig_y, fig_dy, fig_d2y};
ylabels = {'$y(t)$', '$y''(t)$', '$y''''(t)$'};


%% Compute values

[tout, xout, stats, mf, Pf, model] = odeFilter(ode, tspan, y0, opts);
filtereval_n = @(t,n) filtereval(t, tout(1:n), mf(:,:,1:n), Pf(:,:,:,1:n), model);

%% First plot -- predict 4

Ncur = 4;

for d=1:3
  figure(figs{d});
  clf;
  
  % Plot function
  plot(t, funcs{d}(t), 'Color',gra, 'LineWidth', 1);
  hold on;
  plot([t0, t0],[-10, 10], '-k', 'LineWidth', 1);
  plot([T, T],[-10, 10], '-k', 'LineWidth', 1);
  
  % Plot prediction/updates n < Ncur
  for n=1:Ncur-1
    tn = t(tout(n) <= t & t <= tout(n+1));
    [mft_n, Pft_n] = filtereval_n(tn, n);
    
    plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    
    mfb_n = filtereval_n([tout(n), tout(n+1)], n);
    plot(tout(n),   mfb_n(d,:,1), '.', 'Color',dre);
    plot(tout(n+1), mfb_n(d,:,2), 'o', 'Color',dre);
  end
  
  % Plot prediction after Ncur
  tn = t(tout(Ncur-1) <= t);
  [mft_n, Pft_n] = filtereval_n(tn, Ncur-1);
  plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  
  % Plot data
  if d == 1
    plot(tout(1), y0, '.', 'Color', ora);
  elseif d == 2
    plot(tout(1:Ncur-1), model.data(1:Ncur-1), '.', 'Color', ora);
  end
  
  % Labels, etc.
  xlim([min(t), max(t)]);
  if d == 1
    ylim([0, 1.]);
  elseif d == 2
    ylim([0, 1.]);  
  else
    ylim([-1.25, 1.25]);
  end 
  
  ylabel(ylabels{d});
  set(gca, 'YTick', []);
  
  set(gca, 'XTick', linspace(t0, T, N+1));
  if d == 3
  xlabel('$t$');
  set(gca, 'XTickLabel', ...
           {'$t_0$', '', '$t_n$', '', '$t_{N-1}$' ,'$T$'});
  else
    set(gca,'XTickLabel',{});
  end
 
end

%% Second plot -- update 4

Ncur = 4;

for d=1:3
  figure(figs{d});
  clf;
  
  % Plot function
  plot(t, funcs{d}(t), 'Color',gra, 'LineWidth', 1);
  hold on;
    plot([t0, t0],[-10, 10], '-k', 'LineWidth', 1);
  plot([T, T],[-10, 10], '-k', 'LineWidth', 1);
  
  % Plot prediction/updates n < Ncur
  for n=1:Ncur-1
    tn = t(tout(n) <= t & t <= tout(n+1));
    [mft_n, Pft_n] = filtereval_n(tn, n);
    
    plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    
    mfb_n = filtereval_n([tout(n), tout(n+1)], n);
    plot(tout(n),   mfb_n(d,:,1), '.', 'Color',dre);
    plot(tout(n+1), mfb_n(d,:,2), 'o', 'Color',dre);
  end
  
  tn = t(tout(Ncur-1) <= t);
  [mft_n, Pft_n] = filtereval_n(tn, Ncur-1);
  plot(tn, squeeze(mft_n(d,:,:)), '--', 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), '--', 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), '--', 'Color',dre);
  
  % Plot update after Ncur
  tn = t(tout(Ncur) <= t);
  [mft_n, Pft_n] = filtereval_n(tn, Ncur);
  plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  
  plot(tout(Ncur), mf(d,:,Ncur), '.', 'Color', dre);
  
  % Plot data
  if d == 1
    plot(tout(1), y0, '.', 'Color', ora);
  elseif d == 2
    plot(tout(1:Ncur-1), model.data(1:Ncur-1), '.', 'Color', ora);
  end
  
  % Labels, etc.
    xlim([min(t), max(t)]);
  if d == 1
    ylim([0, 1.]);
  elseif d == 2
    ylim([0, 1.]);  
  else
    ylim([-1.25, 1.25]);
  end 
  
  % ylabel(ylabels{d});
  set(gca, 'YTick', []);
  
  set(gca, 'XTick', linspace(t0, T, N+1));
  if d == 3
  xlabel('$t$');
  set(gca, 'XTickLabel', ...
           {'$t_0$', '', '$t_n$', '', '$t_{N-1}$' ,'$T$'});
  else
    set(gca,'XTickLabel',{});
  end
  
end

%% Third plot -- predict 5

Ncur = 5;

for d=1:3
  figure(figs{d});
  clf;
  
  % Plot function
  plot(t, funcs{d}(t), 'Color',gra, 'LineWidth', 1);
  hold on;
    plot([t0, t0],[-10, 10], '-k', 'LineWidth', 1);
  plot([T, T],[-10, 10], '-k', 'LineWidth', 1);
  
  % Plot prediction/updates n < Ncur
  for n=1:Ncur-1
    tn = t(tout(n) <= t & t <= tout(n+1));
    [mft_n, Pft_n] = filtereval_n(tn, n);
    
    plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    
    mfb_n = filtereval_n([tout(n), tout(n+1)], n);
    plot(tout(n),   mfb_n(d,:,1), '.', 'Color',dre);
    plot(tout(n+1), mfb_n(d,:,2), 'o', 'Color',dre);
  end
  
  % Plot prediction after Ncur
  tn = t(tout(Ncur-1) <= t);
  [mft_n, Pft_n] = filtereval_n(tn, Ncur-1);
  plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  
  % Plot data
  if d == 1
    plot(tout(1), y0, '.', 'Color', ora);
  elseif d == 2
    plot(tout(1:Ncur-1), model.data(1:Ncur-1), '.', 'Color', ora);
  end
  
  % Labels, etc.
  xlim([min(t), max(t)]);
  if d == 1
    ylim([0, 1.]);
  elseif d == 2
    ylim([0, 1.]);  
  else
    ylim([-1.25, 1.25]);
  end 
  
  % ylabel(ylabels{d});
  set(gca, 'YTick', []);
  
  set(gca, 'XTick', linspace(t0, T, N+1));
  if d == 3
  xlabel('$t$');
  set(gca, 'XTickLabel', ...
           {'$t_0$', '', '$t_n$', '', '$t_{N-1}$' ,'$T$'});
  else
    set(gca,'XTickLabel',{});
  end
  
end

%% Fourth plot -- update 5

Ncur = 5;

for d=1:3
  figure(figs{d});
  clf;
  
  % Plot function
  plot(t, funcs{d}(t), 'Color',gra, 'LineWidth', 1);
  hold on;
    plot([t0, t0],[-10, 10], '-k', 'LineWidth', 1);
  plot([T, T],[-10, 10], '-k', 'LineWidth', 1);
  
  % Plot prediction/updates n < Ncur
  for n=1:Ncur-1
    tn = t(tout(n) <= t & t <= tout(n+1));
    [mft_n, Pft_n] = filtereval_n(tn, n);
    
    plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
    
    mfb_n = filtereval_n([tout(n), tout(n+1)], n);
    plot(tout(n),   mfb_n(d,:,1), '.', 'Color',dre);
    plot(tout(n+1), mfb_n(d,:,2), 'o', 'Color',dre);
  end
  
  tn = t(tout(Ncur-1) <= t);
  [mft_n, Pft_n] = filtereval_n(tn, Ncur-1);
  plot(tn, squeeze(mft_n(d,:,:)), '--', 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), '--', 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), '--', 'Color',dre);
  
  tn = t(tout(Ncur) <= t);
  [mft_n, Pft_n] = filtereval_n(tn, Ncur);
  plot(tn, squeeze(mft_n(d,:,:)), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) + 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  plot(tn, squeeze(mft_n(d,:,:)) - 2*sqrt(squeeze(Pft_n(d,d,:,:))), 'Color',dre);
  
  plot(tout(Ncur), mf(d,:,Ncur), '.', 'Color', dre);
  
  % Plot data
  if d == 1
    plot(tout(1), y0, '.', 'Color', ora);
  elseif d == 2
    plot(tout(1:Ncur-1), model.data(1:Ncur-1), '.', 'Color', ora);
  end
  
  % Labels, etc.
  xlim([min(t), max(t)]);
  if d == 1
    ylim([0, 1.]);
  elseif d == 2
    ylim([0, 1.]);  
  else
    ylim([-1.25, 1.25]);
  end 
  
  % ylabel(ylabels{d});
  set(gca, 'YTick', []);
  
  set(gca, 'XTick', linspace(t0, T, N+1));
  if d == 3
  xlabel('$t$');
  set(gca, 'XTickLabel', ...
           {'$t_0$', '', '$t_n$', '', '$t_{N-1}$' ,'$T$'});
  else
    set(gca,'XTickLabel',{});
  end
  
end

%% Fifth plot -- smoothed and sampled

[ms, Ps] = odeSmoother(tout, mf, Pf, model);

tf = t(tout(1) <= t);
[mst, Pst] = odeSmoother(tout, mf, Pf, model, num_samples, tf);
[~, ~, St] = odeSmoother(tout, mf, Pf, model, num_samples, tf(1:4:end));

for d=1:3
  figure(figs{d});
  clf;
  
  % Plot function
  plot(t, funcs{d}(t), 'Color',gra, 'LineWidth', 1);
  hold on;
    plot([t0, t0],[-10, 10], '-k', 'LineWidth', 1);
  plot([T, T],[-10, 10], '-k', 'LineWidth', 1);
  
  plot(tf, squeeze(mst(d,:,:)), 'Color',mpg);
  plot(tf, squeeze(mst(d,:,:)) + 2*sqrt(squeeze(Pst(d,d,:,:))), 'Color',mpg);
  plot(tf, squeeze(mst(d,:,:)) - 2*sqrt(squeeze(Pst(d,d,:,:))), 'Color',mpg);

  for k=1:num_samples
    plot(tf(1:4:end), squeeze(St(d,:,:,k))', '-', 'Color',lightmpg);
  end

  % Plot data
  if d == 1
    plot(tout(1), y0, '.', 'Color', ora);
  elseif d == 2
    plot(tout, model.data, '.', 'Color', ora);
  end
  
  % Labels, etc.
   xlim([min(t), max(t)]);
  if d == 1
    ylim([0, 1.]);
  elseif d == 2
    ylim([0, 1.]);  
  else
    ylim([-1.25, 1.25]);
  end 
  
  % ylabel(ylabels{d});
  set(gca, 'YTick', []);
  
  set(gca, 'XTick', linspace(t0, T, N+1));
  if d == 3
  xlabel('$t$');
  set(gca, 'XTickLabel', ...
           {'$t_0$', '', '$t_n$', '', '$t_{N-1}$' ,'$T$'});
  else
    set(gca,'XTickLabel',{});
  end
  
end

