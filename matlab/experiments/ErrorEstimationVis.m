% ErrorEstimationVis.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2017-07-13
% Version: 0.1
% Purpose: Produces histograms of [est. error] / [true error]

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

style = 0; % two variants: plot 5 problems in one plot or plot each problem

num_bins    = 16;
max_c       = 1.6;
bin_centers = linspace(0 + max_c / (2*num_bins), max_c - max_c / (2*num_bins), num_bins);

ts = linspace(0, max_c, 100);
chicdf = 0.5 * (erf(ts./sqrt(2)) - erf(-ts./sqrt(2)));

%% Configuration

tol = 1e-3;
% tol = [1e-3, 1e-6];
% tol = 1e-6;
% tol = [1e-3, 1e-6, 1e-9];
% tol = 1e-9;

odes = hullBenchmark;
ode_sel = 1:numel(odes);

options = struct ();
options.N = 2;

options.alpha  = 0 * 0.001;      % how much should the uncertainty of P be weighted
options.lambda = 0.8;     % learning rate

options.ErrorPerUnitStep = true;
options.eta = [0.1, 5.0]; % min, max in h = eta * h
options.rho = 0.95;        % consider rho * tol for safety
options.nu  = 0.75;        % h = nu * h if step failed
options.k_I = (options.N + 1) / 2;   % weight of the integration controller
% options.k_I = 2 / 2;
options.k_P = 0;   % weight of the proportional controller

options.InitialStep = 1e-7;

F = diag(1:(options.N),1);
L = [zeros(options.N,1); 1/factorial(options.N)];

all_counts = cell(numel(ode_sel), numel(tol));

% for i=1:numel(odes)
for i=1:numel(ode_sel)
  
  info(sprintf('solving problem with idx %02i, %02i/%02i', ...
               ode_sel(i), i,numel(ode_sel)));
  
  ode = odes{ode_sel(i)};
  
  for t=1:numel(tol)
    
    info(sprintf('tolerance lvl %.1e',tol(t)),1);
    
    thisTol = tol(t);
    
    options.ewt = @(y) 1 ./ (0. * abs(y) + thisTol);
    
    [tout, xout, stats, mf, Pf, model] = ...
      odeFilter(ode.fun, ode.tspan, ode.x0, options);
    
    statsHull = evalSolHull(ode, tout, xout, thisTol, options.ErrorPerUnitStep);
    
    % error estimation only works from k-index 3!
    trueErrors = max(abs(statsHull{1}(3:end,:) - xout(3:end,:)),eps);
    
    estErrors = NaN(size(trueErrors));
    for k=3:numel(tout)
       deltaT = tout(k) - tout(k-1);
       A  = sdePhi11Matrix(F, deltaT);
       Q1 = sdePhi12Matrix(F, L, 1, deltaT) * A.';
        
       for d=1:size(xout,2)
           estErrors(k-2,d) = sqrt(model.ssqs(d,k) * Q1(1,1));
       end
    end
    
    relErrors = trueErrors ./ estErrors;
    all_counts{i,t} = hist(relErrors, bin_centers);
    
  end
  
end

% keyboard;

%% Plotting

bin_c_zero = [0, bin_centers];

% prepare data
if style == 1
    
    if size(all_counts,1) ~= 25
        error('this assumes you have solved all problems!');
    end
    
    new_all_counts = cell(5,numel(tol));
    for t=1:numel(tol)
        for c=1:5
            sub_counts = NaN(numel(bin_centers), 5);
            for p=1:5
                this_counts = all_counts{(c-1)*5 + p, t};
                if size(this_counts,1) > 1
                    this_counts = sum(this_counts,2);
                else
                    this_counts = this_counts';
                end
                sub_counts(:,p) = this_counts;
            end
            new_all_counts{c,t} = sub_counts;
        end
    end
    
    all_counts = new_all_counts;
    
else
    
    for t=1:numel(tol)
       for c=1:numel(ode_sel)
           if size(all_counts{c,t},1) == 1
               all_counts{c,t} = all_counts{c,t}'; 
           end
       end
    end
    
end

for t=1:numel(tol)
    
    for c=1:size(all_counts,1)
    
        thisplot = sprintf('style%1i-set%02i-tol%1.0e', style, c, tol(t));
        
        clf;
        hold on;
        
        this_counts = all_counts{c,t};
        
        counts = [this_counts(1,:); this_counts];
        
        col_grad = linspace(0., 1., size(counts,2)).^0.5;
        cols     = bsxfun(@times, blu, col_grad.') + bsxfun(@times, lightblu, (1 - col_grad).');
        
        for d=1:size(counts,2)
            this_col = cols(d,:);
            plot(bin_c_zero, cumsum(counts(:,d)) ./ sum(counts(:,d)), 'color', this_col);
        end
        
        box on;
        xlim([0, max_c]);
        ylim([0, 1]);
        
        plot(ts, chicdf, '--', 'color', mpg);
        plot([1, 1], [0, 1], 'Color', dre);
        
        set(gca, 'YTick', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]);
        set(gca, 'YTickLabel', {});
        set(gca, 'XTick', [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.55]);
        set(gca, 'XTickLabel', {0, '', 0.5, '', 1, '', '', '$>1.5$'});
        
        drawnow;
        keyboard;
   
    end
    
end