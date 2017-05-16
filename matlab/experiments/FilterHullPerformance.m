% FilterHullPerformance.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2016-05-25
% Version: 0.1
% Purpose: Test the filter on the Hull benchmark

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

info('This takes a long time');
info('Computation can be speed up if the parallel processing toolbox is available')
info('In this case: please change the for-loop in odeSeqGT to parfor')

% tol = 1e-3;
tol = [1e-3, 1e-6];
% tol = [1e-3, 1e-6, 1e-9];

odes = hullBenchmark;
ode_sel = 1:numel(odes);

options = struct ();
options.N = 2;

options.alpha  = 0 * 0.001; % how much should the uncertainty of P be weighted
options.lambda = 0.8;       % learning rate

options.ErrorPerUnitStep = true;
options.eta = [0.1, 5.0]; % min, max in h = eta * h
options.rho = 0.95;        % consider rho * tol for safety
options.nu  = 0.75;        % h = nu * h if step failed
options.k_I = (options.N + 1) / 2;   % weight of the integration controller
options.k_P = 0;   % weight of the proportional controller

options.InitialStep = 1e-7;

num_samples = 0;

fcnCalls = NaN(numel(ode_sel), numel(tol));
numSteps = NaN(numel(ode_sel), numel(tol));
percDecv = NaN(numel(ode_sel), numel(tol));
maxError = NaN(numel(ode_sel), numel(tol));

for i=1:numel(ode_sel)
  
  info(sprintf('solving problem with idx %02i, %02i/%02i', ...
               ode_sel(i), i,numel(ode_sel)));
  
  ode = odes{ode_sel(i)};
  
  for t=1:numel(tol)
    
    info(sprintf('tolerance lvl %.1e',tol(t)),1);
    
    thisTol = tol(t);
    
    options.ewt = @(y) 1 ./ (0. * abs(y) + thisTol);
    
    [tout, xout, stats] = ...
      odeFilter(ode.fun, ode.tspan, ode.x0, options);
    
    statsHull = evalSolHull(ode, tout, xout, thisTol, options.ErrorPerUnitStep);
    
    fcnCalls(i,t) = stats.fcnCalls;
    numSteps(i,t) = statsHull{2};
    
    percDecv(i,t) = statsHull{3};
    maxError(i,t) = statsHull{4};
    
  end
  
end

%% Plotting

vcmap = viridis(128);

figFcnCalls = figure();
load('solver/benchmark/data/fcnCalls.mat');
thisData = fcnCalls;

solvers = {'Extrap.','Krogh','Gear','RK4','RK6','RK8','ode23','ode45','ode113','PNM'};

for t=1:numel(tol)
  data_size = size(data);
  subdata = reshape(data(:,:,:,t),[data_size(1), data_size(2) * data_size(3)]);
  subdata = subdata(:,ode_sel);
  
  subplot(1,numel(tol),t);
  hold on;
  imagesc(log10([subdata; thisData(:,t).']));
  
  if t == 1
    set(gca,'YTick',1:10);
    set(gca,'YTickLabels',solvers);
  else
    set(gca,'YTick',[]);
  end
  
  set(gca,'XTick',2.5:5:22.5)
  set(gca,'XTickLabels',{'Single', 'Small','Linear','Orbit','Higher'});
  
  plot([0.5,25.5],[9.5, 9.5],'-w','LineWidth',1);
  for i=1:4
    plot([5*i, 5*i] + .5, [0.5, 10.5],'-w','LineWidth',0.5);
  end
  
  colormap(vcmap);
  colorbar;
  % caxis([1.5, 4]);
  % lcolorbar('$\log_{10}$(#FE)');
  title(sprintf('Tol: %.1e',tol(t)));
  axis ij tight;
end



figNumSteps = figure();
load('solver/benchmark/data/numSteps.mat');
thisData = numSteps;

solvers = {'Extrap.','Krogh','Gear','RK4','RK6','RK8','ode23','ode45','ode113','PNM'};

for t=1:numel(tol)
  data_size = size(data);
  subdata = reshape(data(:,:,:,t),[data_size(1), data_size(2) * data_size(3)]);
  subdata = subdata(:,ode_sel);
  
  subplot(1,numel(tol),t);
  hold on;
  imagesc(log10([subdata; thisData(:,t).']));
  
  if t == 1
    set(gca,'YTick',1:10);
    set(gca,'YTickLabels',solvers);
  else
    set(gca,'YTick',[]);
  end
  
  set(gca,'XTick',2.5:5:22.5)
  set(gca,'XTickLabels',{'Single', 'Small','Linear','Orbit','Higher'});
  
  plot([0.5,25.5],[9.5, 9.5],'-w','LineWidth',1);
  for i=1:4
    plot([5*i, 5*i] + .5, [0.5, 10.5],'-w','LineWidth',0.5);
  end
  
  colormap(vcmap);
  colorbar;
  % caxis([0.5, 3.5]);
  title(sprintf('Tol: %.1e',tol(t)));
  axis ij tight;
end



figPercDecv = figure();
load('solver/benchmark/data/percDecv.mat');
thisData = percDecv;

solvers = {'Extrap.','Krogh','Gear','RK4','RK6','RK8','ode23','ode45','ode113','PNM'};

for t=1:numel(tol)
  data_size = size(data);
  subdata = reshape(data(:,:,:,t),[data_size(1), data_size(2) * data_size(3)]);
  subdata = subdata(:,ode_sel);
  
  subplot(1,numel(tol),t);
  hold on;
  imagesc(log10([subdata; thisData(:,t).']));
  
  if t == 1
    set(gca,'YTick',1:10);
    set(gca,'YTickLabels',solvers);
  else
    set(gca,'YTick',[]);
  end
  
  set(gca,'XTick',2.5:5:22.5)
  set(gca,'XTickLabels',{'Single', 'Small','Linear','Orbit','Higher'});
  
  plot([0.5,25.5],[9.5, 9.5],'-w','LineWidth',1);
  for i=1:4
    plot([5*i, 5*i] + .5, [0.5, 10.5],'-w','LineWidth',0.5);
  end
  
  colormap(vcmap);
  colorbar;
  % caxis([-3,1.25]);
  title(sprintf('Tol: %.1e',tol(t)));
  axis ij tight;
end



figMaxError = figure();
load('solver/benchmark/data/maxError.mat');
thisData = maxError;

solvers = {'Extrap.','Krogh','Gear','RK4','RK6','RK8','ode23','ode45','ode113','PNM'};

for t=1:numel(tol)
  data_size = size(data);
  subdata = reshape(data(:,:,:,t),[data_size(1), data_size(2) * data_size(3)]);
  subdata = subdata(:,ode_sel);
  
  subplot(1,numel(tol),t);
  hold on;
  imagesc(log10([subdata; thisData(:,t).']));
  
  if t == 1
    set(gca,'YTick',1:10);
    set(gca,'YTickLabels',solvers);
  else
    set(gca,'YTick',[]);
  end
  
  set(gca,'XTick',2.5:5:22.5)
  set(gca,'XTickLabels',{'Single', 'Small','Linear','Orbit','Higher'});
  
  plot([0.5,25.5],[9.5, 9.5],'-w','LineWidth',1);
  for i=1:4
    plot([5*i, 5*i] + .5, [0.5, 10.5],'-w','LineWidth',0.5);
  end
  
  colormap(vcmap);
  colorbar;
  % caxis([-1.5,2]);
  title(sprintf('Tol: %.1e',tol(t)));
  axis ij tight;
end
