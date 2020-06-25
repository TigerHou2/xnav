close all
clear;clc

% define orbit
mu = 1;
a = 1;
e = 0.2;
i = pi/4;
omg = 3*pi/2;
w = pi/3;
params = [a e i omg w 0];

% define time offset
t_offset = 3.9; % just has to be a number unrelated to anything else
period = 2*pi * sqrt(a^3/mu);
t_offset = mod(t_offset,period);

% define pulsars
P = [0 1 1;... % pulsar 1
     1 1 1;... % pulsar 2
     3 1 2]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% calculate true position and velocity values
% --- at least 3 measurements for each pulsar
% --- if some pulsars have fewer measurements, fill to end with NaN
f = [ 0  40  80;... % pulsar 1
     10  50  90;... % pulsar 2
     20  60  100]; % pulsar 3
f = f-40;
f = deg2rad(f);
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
t_true = sqrt(a^3/mu)*M;

t = t_true + t_offset;

% perform observations
r = nan(size(f,1),size(f,2),3);
v = nan(size(f,1),size(f,2),3);
obsv = nan(size(f));
pulsar = P;
for i = 1:size(f,1)
    for j = 1:size(f,2)
        params(6) = f(i,j);
        [r(i,j,:),v(i,j,:)] = Get_Orb_Vects(params,mu);
%         noise = randn(1,1,3);
%         noise = noise ./ vecnorm(noise,2,3) * 0.0005;
%         vtemp = v(i,j,:) + noise;
        vtemp = v(i,j,:);
        obsv(i,j) = P(:,i)'*vtemp(:);
    end
end

% visualize range rate
figure(1)
scatter(f(:),obsv(:))

% calculate perfect inputs
OPT = [e,period,t_offset];

% check if perfect input yields perfect output
[optDiff,V] = rrFun_sine(OPT,obsv,pulsar,mu,t)

%% global search

warning('off','all')

res = [15,15,25];
dat = nan(res);
ee = linspace(0,0.9,res(1));
pmin = max(t(:))-min(t(:));
pmax = 10*(max(t(:))-min(t(:)));
pp = linspace(pmin,pmax,res(2));
tt = linspace(0,pmax,res(3));

E = nan(res(1),res(2),res(3));
P = nan(res(1),res(2),res(3));
T = nan(res(1),res(2),res(3));

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    if tt(k) > pp(j)
        % we don't care about redunrant guesses where the time since
        % periapsis is greater than the period
        dat(i,j,k) = Inf;
    else
        fin = [ee(i),pp(j),tt(k)];
        dat(i,j,k) = norm(rrFun_sine(fin,obsv,pulsar,mu,t));
    end
    E(i,j,k) = ee(i);
    P(i,j,k) = pp(j);
    T(i,j,k) = tt(k);
end
end
end

%----- testing -----
% use convolution to find element with best neighbors
% --- version 1: poor man's gaussian blur
%     neighbors = 2;
%     edge = 2*neighbors + 1;
%     filter = ones(edge,edge,edge);
%     dat = convn(dat,filter,'same') + dat*edge^3;
% --- version 2: gaussian blur
    dat = imgaussfilt3(dat);

[~,idx] = min(dat(:));
[i,j,k] = ind2sub(size(dat),idx);
f0 = [ee(i(1)),pp(j(1)),tt(k(1))];

% scatter to see solution space
figure(2)
Emesh = E(:);
Pmesh = P(:);
Tmesh = T(:);
D = dat(:).^0.25 * 15;
cutoff = 40;
Emesh = Emesh(D<cutoff);
Pmesh = Pmesh(D<cutoff);
Tmesh = Tmesh(D<cutoff);
D = D(D<cutoff);
scatter3(OPT(1),OPT(2),OPT(3),24,'magenta','filled','DisplayName','True Soln')
hold on
scatter3( f0(1), f0(2), f0(3),24,'red','DisplayName','Init Guess')
f1 = scatter3(Emesh,Pmesh,Tmesh,18,D,'filled','DisplayName','Objective Function');
f1.MarkerFaceAlpha = 0.6;
hold off
colormap(bone(200))
xlabel('eccentricity')
ylabel('period')
zlabel('time since periapse')
pbaspect([1 1 1])
colorbar
legend('Location','Best')

%% optimization

disp('Searching for solution...')

fun = @(x) rrFun_sine(x,obsv,pulsar,mu,t);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',700 ...
                               ,'Algorithm','levenberg-marquardt' ...
                               ,'StepTolerance',1e-8 ...
                               ,'FunctionTolerance',1e-8);
g_opt = fsolve(fun,f0,options);
soln_opt = g_opt;
soln_opt(3) = mod(g_opt(3),g_opt(2))
close(gcf)

[optDiff,V] = rrFun_sine(soln_opt,obsv,pulsar,mu,t,'debug');

%% plot velocities

figure(3)
t1 = V(:,:,1);
t2 = V(:,:,2);
t3 = V(:,:,3);
scatter3(t1(:),t2(:),t3(:),'DisplayName','Guess Vel')
hold on
scatter3(t1(1),t2(1),t3(1),'Filled','DisplayName','Guess Vel Initial')
t1 = v(:,:,1);
t2 = v(:,:,2);
t3 = v(:,:,3);
scatter3(t1(:),t2(:),t3(:),'DisplayName','True Vel')
scatter3(t1(1),t2(1),t3(1),'Filled','DisplayName','True Vel Initial')
hold off
axis equal
legend('Location','Best')