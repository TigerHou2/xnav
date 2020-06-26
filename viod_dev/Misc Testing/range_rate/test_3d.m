close all
clear;clc

mu = 1;
a = 1;
e = 0.8;
i = pi/5;
omg = 3*pi/2;
w = pi/3;

params = [a e i omg w 0];

% define pulsars
P = [0 1 1;... % pulsar 1
     1 1 1;... % pulsar 2
     3 1 2]';  % pulsar 3
P = P ./ vecnorm(P,2,1);

% calculate true position and velocity values
f = deg2rad([0 10 20 40 60 70 80]);
E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M = E - e*sin(E);
t = sqrt(a^3/mu)*M;
dt = t(2:end) - t(1:end-1);
r = nan(3,length(f));
v = nan(3,length(f));
for j = 1:length(f)
    params(6) = f(j);
    [r(:,j),v(:,j)] = Get_Orb_Vects(params,mu);
end

obsv = nan(1,length(f));
pulsar = nan(3,length(f));
% alternate between three pulsars
idx = 0;
for j = 1:length(f)
    idx = mod(idx,3) + 1;
    obsv(j) = v(:,j)' * P(:,idx);
    pulsar(:,j) = P(:,idx);
end

[~,~,V] = svd(v');
k = V(:,end);

% check if orbit plane normal is in the right direction
% we assume more than half of the measurements are taken less than half an
% orbit apart from its adjacent measurements.
kEst = zeros(3,1);
for i = 1:size(v,2)-1
    kEst = kEst + cross(v(:,i),v(:,i+1));
end
if dot(k,kEst) < 0
    k = -k;
end

ux = cross(v(:,1),k);
ux = ux / norm(ux);
uy = cross(k,ux);
T = [ux';uy';k'];

v

% ----- true hodograph ------
vt = T * v;
[C,R] = fitcircle(vt(1:2,:));
Z = T' * [C;0];
% find theta
uc = Z / norm(Z);
ue = cross([0;0;1],uc);
ue = ue / norm(ue);
k0 = cross(uc,ue);
theta = atan2(norm(cross(k0,k)),dot(k0,k));
if dot(ue,k)>0
    theta = -theta;
end

    % visualize for testing
    figure
    % --- calculate points on hodograph
    angles = linspace(0,2*pi,100);
    cx = R * cos(angles) + C(1);
    cy = R * sin(angles) + C(2);
    circ = T' * [cx;cy;zeros(size(cx))];
    % --- plot hodograph
    plot3(circ(1,:),circ(2,:),circ(3,:),'DisplayName','Hodograph');
    hold on
    % --- plot inertial axes
    quiver3([0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1],'Color','Black')
    % --- plot circle center vector uc
    quiver3(0,0,0,uc(1),uc(2),uc(3),'DisplayName','uc');
    % --- plot eccentricity vector ue
    quiver3(0,0,0,ue(1),ue(2),ue(3),'DisplayName','ue')
    % --- plot hodograph normal at theta = 0
    quiver3(0,0,0,k0(1),k0(2),k0(3),'DisplayName','k0')
    % --- plot true hodograph normal
    quiver3(0,0,0,k(1),k(2),k(3),'DisplayName','k')
    % --- plot pulsar vector
    quiver3(0,0,0,P(1,1),P(2,1),P(3,1),'-.','LineWidth',1,'DisplayName','P1')
    quiver3(0,0,0,P(1,2),P(2,2),P(3,2),'--','LineWidth',1,'DisplayName','P2')
    quiver3(0,0,0,P(1,3),P(2,3),P(3,3),':' ,'LineWidth',1,'DisplayName','P3')
    % --- plot measurement points
    scatter3(v(1,:),v(2,:),v(3,:))
    hold off
    xlabel('x');ylabel('y');zlabel('z');legend
    axis equal
    view([1 1 1])

% test if objective function returns 0 given perfect input
hodo_true = [Z',theta,R];
% [err,val] = rrFun([Z',theta,R],v_obsv,mu,dt)
[err,vel] = rrFun_ta([Z',theta,R],obsv,pulsar,mu,dt);

disp(['Max range-rate error: ' num2str(max(err))])
disp(['Max velocity error:   ' num2str(max(max(abs(v-vel))))])

%% global search

res = [8,8,8,8,8];
dat = nan(res);
lb = mean(abs(obsv));
ub = 3*mean(abs(obsv));
xx = linspace(-lb,lb,res(1)); % hodograph center x pos
yy = linspace(-lb,lb,res(2)); % hodograph center y pos
zz = linspace(-lb,lb,res(3)); % hodograph center z pos
tt = linspace(-pi,pi,res(4)); % hodograph rotation about uc
rr = linspace( lb,ub,res(5)); % hodograph radius
idx = 0;

for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
for m = 1:res(4)
for n = 1:res(5)
    idx = idx + 1;
    fin = [xx(i),yy(j),zz(k),tt(m),rr(n)];
    dat(i,j,k,m,n) = norm(rrFun_ta(fin,obsv,pulsar,mu,dt));
end
end
end
end
end

% use convolution to find element with best neighbors
% --- version 1: poor man's gaussian blur
    neighbors = 2;
    edge = 2*neighbors + 1;
    filter = ones(edge,edge,edge,edge,edge);
    dat = (convn(dat,filter,'same') + dat*edge^5) / edge^5 / 2;

[~,idx] = min(dat(:));
[i,j,k,m,n] = ind2sub(size(dat),idx);
f0 = [xx(i(1)),yy(j(1)),zz(k(1)),tt(m(1)),rr(n(1))];

%% optimization

disp('Searching for solution...')

% f0 = hodo_true .* [1.02 1.02 1.02 1.02 1.02];

f = @(x) rrFun_ta(x,obsv,pulsar,mu,dt);
options = optimoptions('fsolve','Display','iter' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',700 ...
                               ,'Algorithm','levenberg-marquardt' ...
                               ,'StepTolerance',1e-8 ...
                               ,'FunctionTolerance',1e-8);
hodo = fsolve(f,f0,options)

[v_err,vel] = rrFun_ta(hodo,obsv,pulsar,mu,dt,1);

v_diff = (v-vel)./vecnorm(v);
if max(abs(v_diff(:))) > 1e-3
    warning('Velocities do not converge!')
end


%% test solution space

f = @(x) rrFun_ta(x,obsv,pulsar,mu,dt);
options = optimoptions('fsolve','Display','none' ...
                               ,'MaxFunctionEvaluations',3000 ...
                               ,'MaxIterations',500 ...
                               ,'Algorithm','levenberg-marquardt' ...
                               ,'StepTolerance',1e-9);

guess_valid = [];
guess_inval = [];
for i = 1:300
    
    pert = (rand(1,5)-0.5);
    % ----- test section -----
% 	pert(end) = abs(pert(end));
    % ----- end test -----
    f0 = hodo_true + hodo_true .* pert;
    hodo = fsolve(f,f0,options);
    if max(abs(hodo-hodo_true)) < 1e-3
        guess_valid(end+1,:) = pert;
    else
        guess_inval(end+1,:) = pert;
    end
    
end

figure;
scatter(guess_inval(:,4),guess_inval(:,5),20);
hold on
scatter(guess_valid(:,4),guess_valid(:,5),20);
hold off
axis equal
xlabel('x')
ylabel('y')

disp(' ')
disp('Mean Disturbance')
disp(['Inval Case:' mat2str(mean(guess_inval),3)])
disp(['Valid Case:' mat2str(mean(guess_valid),3)])
disp(' ')
disp('Standard Deviation')
disp(['Inval Case:' mat2str(std(guess_inval),3)])
disp(['Valid Case:' mat2str(std(guess_valid),3)])
disp(' ')
disp('Correlation Coefficients:')
disp(corrcoef(guess_valid))