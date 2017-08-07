% Tutorial.m
% Author: Michael Schober (mschober@tue.mpg.de)
% Date: 2017-05-15
% Version: 0.1
% Purpose: Demonstrate basic usage of odeFilter and compare it to ode45

%% Setup

clear;
close all;

setup;


%% Configuration

% define the right-hand-side of the ODE.
% We'll use an anonymous function, but you can also put it in an extra file
% f : [t,y] -> y'
A = -1;
f =  @(t,y) A*y;

% next, define the interval you'll require the solution on
t0 = 0.;
T  = 1.;

tspan = [t0, T];

% finally, define y0 = y(t_0)
y0 = 1;


%% Test the problem with ode45, the standard Matlab solver
% input: function handle, 2-vector with start and end, n-vector of initial
% values
% output: vector with [t0, t1, t2, ..., tN = T], vector with [y0 = y(t0),
% y1 approx y(t1), y2 approx y(t2), ..., yN]
[tout, yout] = ode45(f, tspan, y0);

% plotting
figMain = figure();
plot(tout,yout,'.-'); % we mark the discrete values, but draw a line between
title('y'' = f(t,y) on [t0, T], y(t0) = y0, solved with ode45');
xlabel('t');
ylabel('y');
keyboard(); % introduces a break here


%% Test the problem with odeFilter
% first: same input and output
[toutFilter, youtFilter] = odeFilter(f, tspan, y0);

% plotting
figure(figMain); clf; % draw in the same figure, but remove lines
plot(toutFilter,youtFilter, '.-');
title('y'' = f(t,y) on [t0, T], y(t0) = y0, solved with odeFilter');
xlabel('t');
ylabel('y');
keyboard();


%% Test with fourth-order Runge-Kutta and fixed step-size
% define step-size
h = 1e-3;
% h = 1e-2;

% Matlab does not have a fixed step-size solver.
% We'll call naiveERK instead, a simple implementation of fixed step-size
% explicit Runge-Kutta methods. This requires a Butcher-tableau object.

btClassic = RKMethods.ClassicRK4; % classic fourth-order method
% the only fourth-order available as filter
btFilter  = RKMethods.FourthOrder(1/3, 1/2);

[toutFixed, youtFixed] = naiveERK(f ,tspan, y0, btClassic, h);

% plotting
figure(figMain); clf; % draw in the same figure, but remove lines
plot(toutFixed,youtFixed, '.-');
title('y'' = f(t,y) on [t0, T], y(t0) = y0, solved with fixed step-size RK');
xlabel('t');
ylabel('y');
keyboard();


%% Test with fourth-order filter and fixed step size

% This requires some extra options. Don't worry about it on the first read.
options = struct (); % a struct is the Matlab version of a dictionary
options.InitialStep = h;
options.N = 4; % fourth-order method
% h_new = h_old * max(eta(1), min(eta(2), proposed_factor))
% i.e.: eta(1)*h_old <= h_new <= eta(2)*h_old
options.eta = [1, 1];
% accept every step, cause we don't measure any error (ewt = error weights)
options.ewt = @(y) zeros(size(y));
options.lambda = 0.; % don't adjust ssq: important to get true multistep methods

% output: first two as above, stats = struct containing info about run
%    mf = marginal posterior mean of Kalman filter at toutFF
%    Pf = marginal posterior covariance of Kalman filter at toutFF
%    model = more Kalman filter data required to compute smoothing dist.
[toutFF, youtFF, stats, mf, Pf, model] = ...
    odeFilter(f, tspan, y0, options);

% plotting
figure(figMain); clf; % draw in the same figure, but remove lines
plot(toutFF,youtFF, '.-');
title('y'' = f(t,y) on [t0, T], y(t0) = y0, solved with fixed step-size filter');
xlabel('t');
ylabel('y');
keyboard();

%% Demonstrate extra functionality: show full posterior distribution

S = 5; % number of samples to draw

% call the RTS smoother that computes the global posterior
% input: tout, mf, Pf, model as computed with odeFilter, number of samples
% (optional)
% output: ms = marginal posterior mean at toutFF
%         Ps = marginal posterior covariance at toutFF
% (optional) samples = samples that have been drawn from the posterior dist
[ms, Ps, samples] = odeSmoother(toutFF, mf, Pf, model, S);

% it's a one-dimensional problem, therefore, the filter and smoother return
% results with one unnecessary singelton dimension
ms = squeeze(ms);
Ps = squeeze(Ps);

% plotting
figure(figMain); % clf; % this time, we keep the old plot, but draw additional stuff
hold on;
plot(toutFF,ms(1,:), '.-'); % function value is at first dimension
% draw plus/minus two standard deviations
plot(toutFF,ms(1,:) + 2*sqrt(squeeze(Ps(1,1,:)))', '--r');
plot(toutFF,ms(1,:) - 2*sqrt(squeeze(Ps(1,1,:)))', '--r');
plot(toutFF,squeeze(samples(1,1,:,:)));
title('y'' = f(t,y) on [t0, T], y(t0) = y0, solved with fixed step-size filter');
xlabel('t');
ylabel('y');
keyboard();

% show the derivatives
for d=2:size(ms,1)
    figure(figMain);
    clf; hold on;
    plot(toutFF,ms(d,:), '.-'); % function value is at first dimension
    % draw plus/minus two standard deviations
    plot(toutFF,ms(d,:) + 2*sqrt(squeeze(Ps(d,d,:)))', '--r');
    plot(toutFF,ms(d,:) - 2*sqrt(squeeze(Ps(d,d,:)))', '--r');
    plot(toutFF,squeeze(samples(d,1,:,:)));
    title([num2str(d), '-th derivative of y on [t0, T]']);
    xlabel('t');
    ylabel('y');
    keyboard();
end