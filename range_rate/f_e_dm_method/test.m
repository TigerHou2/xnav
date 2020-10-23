close all
clear;clc

addpath('../fcns_circ')
addpath('../fcns_orb')

% define orbit
mu = 1;
a = 1;
i = deg2rad(23);
omg = deg2rad(10);
w = 0;

% define true solution
e = 0.95;
f0 = deg2rad(13);
dM = deg2rad(59);

OPT = [f0,e,dM];

params = [a e i omg w 0];

% define pulsars
P = [0 -1 1;... % pulsar 1
     1 0 1;... % pulsar 2
     -1 0 3]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% calculate mean and true anomalies of measurements
E0 = 2*atan(sqrt((1-e)/(1+e))*tan(f0/2));
M0 = E0 - e*sin(E0);
M = M0 + linspace(0,dM,length(P(:)));
M = reshape(M,size(P'));
E = kepler(M,e);
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

% calculate time of measurements
t_true = (M-M0)/2/pi*sqrt(a^3/mu);

% set random epoch
t_offset = 3.5; % just has to be a number unrelated to anything else
period = 2*pi * sqrt(a^3/mu);
t_offset = mod(t_offset,period);
t_meas = t_true + t_offset;

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
plot(f(1,:),obsv(1,:),'ro')
hold on
plot(f(2,:),obsv(2,:),'go')
plot(f(3,:),obsv(3,:),'bo')
hold off

% visualize pulsars
figure(100)
quiver3(0,0,0,P(1,1),P(2,1),P(3,1),'DisplayName','Pulsar 1')
hold on
quiver3(0,0,0,P(1,2),P(2,2),P(3,2),'DisplayName','Pulsar 2')
quiver3(0,0,0,P(1,3),P(2,3),P(3,3),'DisplayName','Pulsar 3')
hold off
legend('Location','Best')

% check if perfect input yields perfect output
ref_error = guess(OPT,obsv,pulsar,mu,t_meas,'debug');
ref_error

%% global search

tic

disp('Searching for initial guess...')

warning('off','all')

res = [50,50,150];
dat = nan(res);
range_f0 = linspace(0,2*pi,res(1)+1);
range_f0(end) = [];
range_e = linspace(0,1,res(2)+1);
range_e(end) = [];
range_dM = linspace(0,2*pi,res(3)+1);
range_dM(1) = [];

F = nan(res(1),res(2),res(3));
E = nan(res(1),res(2),res(3));
M = nan(res(1),res(2),res(3));

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    fin = [range_f0(i),range_e(j),range_dM(k)];
    out = guess(fin,obsv,pulsar,mu,t_meas);
    dat(i,j,k) = norm(out(:));
    F(i,j,k) = range_f0(i);
    E(i,j,k) = range_e(j);
    M(i,j,k) = range_dM(k);
end
end
end

% mins = islocalmin(dat);
% [~,idx] = min(dat(:)./(mins(:)+0.01));
[~,idx] = min(dat(:));
[i,j,k] = ind2sub(size(dat),idx);
fun0 = [range_f0(i(1)),range_e(j(1)),range_dM(k(1))];


%% scatter global research results

% scatter to see solution space
figure(2)
Fmesh = F(:);
Emesh = E(:);
Mmesh = M(:);
D = dat(:);
lb = quantile(D,0);
ub = quantile(D,0.001);
Fmesh = Fmesh(D>lb&D<ub);
Emesh = Emesh(D>lb&D<ub);
Mmesh = Mmesh(D>lb&D<ub);
D = D(D>lb&D<ub);
mins = islocalmin(dat);
Floc = F(mins);
Eloc = E(mins);
Mloc = M(mins);
hold on
scatter3(OPT(1),OPT(2),OPT(3),24,'cyan',...
        'DisplayName','True Soln',...
        'LineWidth',1.5)
scatter3( fun0(1), fun0(2), fun0(3),24,'green',...
        'DisplayName','Init Guess',...
        'LineWidth',1.5)
fig = scatter3(Fmesh,Emesh,Mmesh,18,D,'filled',...
        'DisplayName','Objective Function');
fig.MarkerFaceAlpha = 0.6;
% scatter3(Floc(:),Eloc(:),Mloc(:),10,'red','filled',...
%         'DisplayName','Local Minima')
hold off

colormap(bone(200))
xlabel('f0')
ylabel('e')
zlabel('dM')
pbaspect([1 1 1])
view([1 1 1])
colorbar
legend('Location','Best')

fun1 = fun0;

%% find the sine wave

disp('Using fsolve...')

fun = @(x) guess(x,obsv,pulsar,mu,t_meas);
options = optimoptions('fsolve','Display','iter' ...
                               ,'MaxFunctionEvaluations',16000 ...
                               ,'StepTolerance', 1e-9 ...
                               ,'FunctionTolerance', 1e-8 ...
                               ,'MaxIterations',4000);
g_opt = fsolve(fun,fun1,options);

soln_opt = g_opt;
soln_opt(1) = mod(g_opt(1),2*pi);
soln_opt(3) = mod(g_opt(3),2*pi);

error = guess(soln_opt,obsv,pulsar,mu,t_meas,'debug');
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