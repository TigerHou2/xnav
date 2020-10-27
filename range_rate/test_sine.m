close all
clear;clc

% define orbit
mu = 1;
a = 1;
e = 0.9;
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
f = f*2;
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
[optDiff,V] = rrFun_sine(OPT,obsv,pulsar,mu,t,'debug');
optDiff


%% global search

tic

disp('Searching for initial guess...')

warning('off','all')

res = [10,100,100];
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
    out = rrFun_sine(fin,obsv,pulsar,mu,t);
    dat(i,j,k) = norm(out(:));
    E(i,j,k) = ee(i);
    P(i,j,k) = pp(j);
    T(i,j,k) = tt(k);
end
end
end

% use convolution to find element with best neighbors
% --- version 1: poor man's gaussian blur
%     neighbors = 2;
%     edge = 2*neighbors + 1;
%     filter = ones(edge,edge);
%     dat = (convn(dat,filter,'same') + dat*edge^2) / edge^2;
% --- version 2: gaussian blur
%     dat = imgaussfilt3(dat);

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
lb = quantile(D,0.00);
ub = quantile(D,0.005);
Emesh = Emesh(D>lb&D<ub);
Pmesh = Pmesh(D>lb&D<ub);
Tmesh = Tmesh(D>lb&D<ub);
D = D(D>lb&D<ub);
scatter3(OPT(1),OPT(2),OPT(3),24,'magenta','filled','DisplayName','True Soln')
hold on
scatter3( f0(1), f0(2), f0(3),24,'red','DisplayName','Init Guess')
fig = scatter3(Emesh,Pmesh,Tmesh,18,D,'filled','DisplayName','Objective Function');
fig.MarkerFaceAlpha = 0.6;
hold off
colormap(bone(200))
xlabel('eccentricity')
ylabel('period')
zlabel('time since periapse')
pbaspect([1 1 1])
colorbar
legend('Location','Best')

f1 = f0;

%% testing section

% f2 = f0;
% f2(end) = OPT(end);
% fJump = nan(1000,3);
% for i = 1:1000
%     [fVal,~,f2] = rrFun_sine(f2,obsv,pulsar,mu,t);
%     disp(num2str(norm(fVal)))
%     disp(mat2str(f2))
%     disp(' ')
%     fJump(i,:) = f2;
% end
% 
% scatter3(fJump(:,1),fJump(:,2),fJump(:,3))

%% optimization

disp('Searching for solution...')

% fun = @(x) norm(rrFun_sine(x,obsv,pulsar,mu,t));
% options = optimoptions('fmincon','Display','iter' ...
%                                 ,'MaxFunctionEvaluations',2000 ...
%                                 ,'MaxIterations',300);
% f1 = fmincon(fun,f0,-eye(3),[0;0;0],[],[],[],[],[],options);

% fun = @(x) norm(rrFun_sine(x,obsv,pulsar,mu,t));
% f1 = fminsearch(fun,f0);

% fun = @(x) norm(rrFun_sine(x,obsv,pulsar,mu,t));
fun = @(x) rrFun_sine(x,obsv,pulsar,mu,t);
options = optimoptions('fsolve','Display','iter' ...
                               ,'MaxFunctionEvaluations',3000 ...
                               ,'StepTolerance', 1e-8 ...
                               ,'MaxIterations',600);
g_opt = fsolve(fun,f1,options);

soln_opt = g_opt;
soln_opt(3) = mod(g_opt(3),g_opt(2));

[optDiff,V,optOut] = rrFun_sine(soln_opt,obsv,pulsar,mu,t,'debug');

disp(['Input:  ' mat2str(soln_opt,4)])
disp(['Output: ' mat2str(optOut,4)])

toc

warning('on','all')

vDiff = V-v;
if max(vDiff(:)) > 1e-3
    warning('Solution does not converge!')
end

%% refine solution
% if output is sufficiently close to true solution, 
% we should be able to iterate a few more times to get the accuracy higher

f3 = soln_opt;
for i = 1:20
    [~,~,f3] = rrFun_sine(soln_opt,obsv,pulsar,mu,t);
end
f3

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