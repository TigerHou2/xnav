close all
clear;clc

% define orbit
mu = 1;
a = 1;
e = 0.7;
i = 0;
omg = 0;
w = 0;
params = [a e i omg w 0];

% define time offset
t_offset = 3.5; % just has to be a number unrelated to anything else
period = 2*pi * sqrt(a^3/mu);
t_offset = mod(t_offset,period);

% define pulsars
P = [0 -1 1;... % pulsar 1
     1 0 1;... % pulsar 2
     -1 0 3]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% calculate true position and velocity values
% --- at least 3 measurements for each pulsar
% --- if some pulsars have fewer measurements, fill to end with NaN
f = [ 0  20  40;... % pulsar 1
      5  25  45;... % pulsar 2
     10  30  50]; % pulsar 3
f = f + 0;
f = deg2rad(f);
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
t_true = sqrt(a^3/mu)*M;

t = t_true + t_offset;
M_offset = t_offset / sqrt(a^3/mu);

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
plot(f(1,:),obsv(1,:),'r-o')
hold on
plot(f(2,:),obsv(2,:),'g-o')
plot(f(3,:),obsv(3,:),'b-o')
hold off

% visualize pulsars
figure(100)
quiver3(0,0,0,P(1,1),P(2,1),P(3,1),'DisplayName','Pulsar 1')
hold on
quiver3(0,0,0,P(1,2),P(2,2),P(3,2),'DisplayName','Pulsar 2')
quiver3(0,0,0,P(1,3),P(2,3),P(3,3),'DisplayName','Pulsar 3')
hold off
legend('Location','Best')

% calculate perfect inputs
OPT = [e,period,M_offset];

% check if perfect input yields perfect output
ref_error = sine_finder(OPT,obsv,pulsar,mu,t,'debug');
ref_error

%% global search

tic

disp('Searching for initial guess...')

warning('off','all')

res = [75,75,50];
dat = nan(res);
ee = linspace(0,0.99,res(1));
pmin = max(t(:))-min(t(:));
% pmax = 2*pi*mu/min(abs(obsv(:)))^3;
pmax = 10;
% pmin = 0.1 * pmax;
pp = linspace(pmin,pmax,res(2));
tt = linspace(0,2*pi,res(3));

% ee = linspace(0,0.99,res(1));
% pp = linspace(period*0.3,period*2,res(2));
% tt = linspace(0,2*pi,res(3));

E = nan(res(1),res(2),res(3));
P = nan(res(1),res(2),res(3));
T = nan(res(1),res(2),res(3));

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    fin = [ee(i),pp(j),tt(k)];
    out = sine_finder(fin,obsv,pulsar,mu,t);
    dat(i,j,k) = norm(out(:));
    E(i,j,k) = ee(i);
    P(i,j,k) = pp(j);
    T(i,j,k) = tt(k);
end
end
end

% mins = islocalmin(dat);
% [~,idx] = min(dat(:)./(mins(:)+0.01));
[~,idx] = min(dat(:));
[i,j,k] = ind2sub(size(dat),idx);
f0 = [ee(i(1)),pp(j(1)),tt(k(1))];


%% scatter global research results

% scatter to see solution space
figure(2)
Emesh = E(:);
Pmesh = P(:);
Tmesh = T(:);
D = dat(:);
lb = quantile(D,0);
ub = quantile(D,0.005);
Emesh = Emesh(D>lb&D<ub);
Pmesh = Pmesh(D>lb&D<ub);
Tmesh = Tmesh(D>lb&D<ub);
D = D(D>lb&D<ub);
mins = islocalmin(dat);
Eloc = E(mins);
Ploc = P(mins);
Tloc = T(mins);
hold on
scatter3(OPT(1),OPT(2),OPT(3),24,'cyan',...
        'DisplayName','True Soln',...
        'LineWidth',1.5)
scatter3( f0(1), f0(2), f0(3),24,'green',...
        'DisplayName','Init Guess',...
        'LineWidth',1.5)
fig = scatter3(Emesh,Pmesh,Tmesh,18,D,'filled',...
        'DisplayName','Objective Function');
fig.MarkerFaceAlpha = 0.6;
scatter3(Eloc(:),Ploc(:),Tloc(:),10,'red','filled',...
        'DisplayName','Local Minima')
hold off

colormap(bone(200))
xlabel('eccentricity')
ylabel('period')
zlabel('time since periapse')
pbaspect([1 1 1])
view([1 1 1])
colorbar
legend('Location','Best')

f1 = f0;

%% find the sine wave

disp('Using fsolve...')

fun = @(x) sine_finder(x,obsv,pulsar,mu,t);
options = optimoptions('fsolve','Display','iter' ...
                               ,'MaxFunctionEvaluations',10000 ...
                               ,'StepTolerance', 1e-9 ...
                               ,'MaxIterations',2500);
g_opt = fsolve(fun,f1,options);

soln_opt = g_opt;
soln_opt(3) = mod(g_opt(3),2*pi);

error = sine_finder(soln_opt,obsv,pulsar,mu,t,'debug');
error

%% function definitions
function Evect = kepler(Mvect,e)

Evect = nan(size(Mvect));

for i = 1:length(Mvect(:))
    
M = Mvect(i);

if M == 0
    Evect(i) = 0;
    continue
end

d = 1;
n = 5;

if e < 1
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(E-e*sin(E)-M) / ...
            (1-e*cos(E) + sign(1-e*cos(E)) ...
                        * sqrt(abs((n-1)^2*(1-e*cos(E))^2 ...
                                   -n*(n-1)*(E-e*sin(E)-M)*(e*sin(E)))));
        d = abs(E-e*sin(E)-M);
    end
else
    E = mod(M+e/2,2*pi);
    while d > 1e-9
        E = E - n*(e*sinh(E)-E-M) / ...
            (e*cosh(E)-1 + sign(e*cosh(E)-1) ...
                         * sqrt(abs((n-1)^2*(e*cosh(E)-1)^2 ...
                                    -n*(n-1)*(e*sinh(E)-E-M)*(e*sinh(E)))));
        d = abs(e*sinh(E)-E-M);
    end
end

Evect(i) = E;

end

end