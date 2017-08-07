% UQVisualization.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2017-07-12
% Version: 0.1
% Purpose: produces figures to visualize various aspects of the uncertainty

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


plot_tubes   = true;

%% Configuration

problem = 2;

num_samples = 5;

% for computations

switch problem
        
    case 1
        
        ode = BRUSWOD;
        ode.tspan = [0, 10];      
        Ds = [1, 2];
        tols  = [1e-1, 1e-1];
        steps = 120;
        gamma = [1e0, 1e0];
        
        xlimits = [0, 5];
        ylimits = xlimits;
        
    case 2
        
        T = 6.6632868593231301896996820305;
        A = 2.00861986087484313650940188;
        
        ode = vanDerPolODE(1, [A; 0], [0, T]);
        Ds  = [1, 2];
        tols = [1e-1, 1e-1];
        steps = 40;
        gamma = [1e0, 1e0];
        
        ode.name = 'VDPOL';
        
        highQSol = odeGroundTruth(ode);
        ode.sol  = @(t) deval(highQSol, t).';
        
        xlimits = [-3, 3];
        ylimits = xlimits;
        
end

% _THETAC_oordinate_S_
thetacs = linspace(0,2*pi,50);
tcs = linspace(ode.tspan(1), ode.tspan(end), 500);
refsol = ode.sol(tcs);

xlabstr = '$y_1(t)$';
ylabstr = '$y_2(t)$';

%% plot solution only

clf; hold on; box on;

plot(refsol(:,Ds(1)), refsol(:,Ds(2)),'-k','LineWidth',1);

xlim(xlimits);
ylim(ylimits);
xlabel(xlabstr);
ylabel(ylabstr);
set(gca, 'XTick', [-1, 0, 1]);
set(gca, 'YTick', [-1, 0, 1]);

drawnow;

%% solve with automatic step size adaptation

opts = struct ();
opts.RelTol = tols(1);
opts.AbsTol = tols(2);

opts.N = 2;
opts.gamma = gamma(1);

ode45sol = ode45(ode.fun, ode.tspan, ode.x0, opts); yout45 = deval(ode45sol, tcs);
[tout, yout, stats, mf, Pf, model] = odeFilter(ode.fun, ode.tspan, ode.x0, opts);
[ms, Ps, samples] = odeSmoother(tout, mf, Pf, model, num_samples, unique([tcs, tout']));

clf; hold on; box on;

plot(refsol(:,Ds(1)), refsol(:,Ds(2)),'-k','LineWidth',1);

plot(yout(:,Ds(1)), yout(:,Ds(2)), '.', 'color', dre, 'MarkerSize', 5);
plot(squeeze(ms(1,Ds(1),:)), squeeze(ms(1,Ds(2),:)), '-', 'Color', dre);
for s=1:num_samples
   plot(squeeze(samples(1,Ds(1),:,s)), squeeze(samples(1,Ds(2),:,s)), ...
        '--', 'Color', dre); 
end


if plot_tubes
    skips = 5;

    for k=3:skips:size(ms,3)
        Rcs = 2 * sqrt(diag([Ps(1,1,Ds(1),k), Ps(1,1,Ds(2),k)]));
        Ccs = bsxfun(@plus, ms(1,Ds,k)', Rcs' * [cos(thetacs); sin(thetacs)]);
        plot(Ccs(1,:), Ccs(2,:),'-', 'Color', mpg);
    end
end

xlim(xlimits);
ylim(ylimits);
xlabel(xlabstr);
ylabel(ylabstr);
set(gca, 'XTick', [-1, 0, 1]);
set(gca, 'YTick', [-1, 0, 1]);

drawnow;
keyboard();

%% solve with fixed step size

opts = struct ();

opts.N = 2;
opts.gamma = gamma(2);

opts.ewt = @(y) zeros(size(y));
opts.eta = [1, 1];
opts.InitialStep = (ode.tspan(end) - ode.tspan(1)) / steps;

[tout, yout, stats, mf, Pf, model] = odeFilter(ode.fun, ode.tspan, ode.x0, opts);
[ms, Ps, samples] = odeSmoother(tout, mf, Pf, model, num_samples, unique([tcs, tout']));

clf; hold on; box on;

plot(refsol(:,Ds(1)), refsol(:,Ds(2)),'-k','LineWidth',1);

plot(yout(:,Ds(1)), yout(:,Ds(2)), 'x', 'color', dre, 'MarkerSize', 5);
plot(squeeze(ms(1,Ds(1),:)), squeeze(ms(1,Ds(2),:)), '-', 'Color', dre);
for s=1:num_samples
    plot(squeeze(samples(1,Ds(1),:,s)), squeeze(samples(1,Ds(2),:,s)), ...
         '--', 'Color', dre); 
end
 
if plot_tubes
    skips = 5;
    for k=3:skips:size(ms,3)
        Rcs = 2 * sqrt(diag([Ps(1,1,Ds(1),k), Ps(1,1,Ds(2),k)]));
        Ccs = bsxfun(@plus, ms(1,Ds,k)', Rcs' * [cos(thetacs); sin(thetacs)]);
        plot(Ccs(1,:), Ccs(2,:),'-', 'Color', mpg);
    end
end

xlim(xlimits);
ylim(ylimits);
xlabel(xlabstr);
ylabel(ylabstr);
set(gca, 'XTick', [-1, 0, 1]);
set(gca, 'YTick', [-1, 0, 1]);

drawnow;
keyboard();

%% use Chkrebtii solver

num_samples_chkr = 20;
   
clf; hold on; box on;

plot(refsol(:,Ds(1)), refsol(:,Ds(2)),'-k','LineWidth',1);

[~, ~, samples] = odeSmoother(tout, mf, Pf, model, num_samples_chkr);

% samples  = sampleFilter(tout, tout(1:end-1), mf(:,:,1:end-1), Pf(:,:,:,1:end-1), model, num_samples_chkr);

for s=1:num_samples_chkr
    plot(squeeze(samples(1,Ds(1),:,s)), squeeze(samples(1,Ds(2),:,s)), ...
         '--', 'Color', dre); 
end

opts.sampled = true;
opts.lambda  = 0.;

for n=1:num_samples_chkr
   [tout_c, xout_c, stats_c, mf_c, Pf_c, model_c] = ...
       odeFilter(ode.fun, ode.tspan, ode.x0, opts);
   
   plot(xout_c(:,Ds(1)), xout_c(:,Ds(2)),'--','Color', mpg);
end

xlim(xlimits);
ylim(ylimits);
xlabel(xlabstr);
ylabel(ylabstr);
set(gca, 'XTick', [-1, 0, 1]);
set(gca, 'YTick', [-1, 0, 1]);

drawnow();
keyboard();
