close all
clear;clc

mu = 1;
a = 1;
e = 0.3;
i = 0;
omg = 0;
w = 0;

params = [a e i omg w 0];

up = [1;0;0]; % pulsar direction
up = up / norm(up);
k = [0;0;1]; % orbit normal

% calculate true values
f = deg2rad([20 55 72 160]);
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

v_obsv = v'*up;

% get the actual hodograph
[Z,R] = fitcircle(v(1:2,:));

f0 = [Z' R/100];
fun(f0,up,v_obsv,mu,dt,k)

res = 300;
dat = nan(res^2,2);
xx = linspace(-0.5,0.5,res);
yy = linspace(-0.5,0.5,res);
[X,Y] = meshgrid(xx,yy);
zz = nan(res,res);
idx = 0;

for i = 1:res
    for j = 1:res
        idx = idx + 1;
        dat(idx,1) = xx(i);
        dat(idx,2) = yy(j);
        zz(i,j) = norm(fun([xx(i),yy(j),R/100],up,v_obsv,mu,dt,k));
    end
end

contour(X,Y,zz)
axis equal

[row,col] = find(zz==min(zz(:)));

%%
f0 = [0.018 0.09 0.02];

f = @(x) fun(x,up,v_obsv,mu,dt,k);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',1000);
hodo = fsolve(f,f0,options)

fun(hodo,up,v_obsv,mu,dt,k)

%% function definitions
function f_diff = fun(C,up,v,mu,dt,k)

    % we will assume 2D for now

    p = C(1); % x-offset of hodograph
    q = C(2); % y-offset of hodograph
    R = C(3)*100; % radius of hodograph
    
    e = sqrt(p^2 + q^2) / R; % eccentricity
    uu = cross(k,up); % vector orthogonal to pulsar direction
    
    % solve intersection on the circle
    
    % --- calculate point on each solution set
    P1 = up * v';
    P2 = P1 + uu;
    
    I1 = nan(2,length(v));
    I2 = nan(2,length(v));
    for j = 1:length(v)
        dx = P2(1,j) - P1(1,j);
        dy = P2(2,j) - P1(2,j);
        dr = sqrt(dx^2+dy^2);
        D = (P1(1,j)-p)*(P2(2,j)-q) - (P2(1,j)-p)*(P1(2,j)-q);
        sgn = sign(dy);
        if sgn == 0
            sgn = 1;
        end
        delta = R^2*dr^2-D^2;
        if delta < 0
            f_diff = delta*ones(size(dt));
            return
        end
        I1(1,j) = ( D*dy + sgn*dx*sqrt(delta) ) / dr^2 + p;
        I1(2,j) = (-D*dx + sgn*dy*sqrt(delta) ) / dr^2 + q;
        I2(1,j) = ( D*dy - sgn*dx*sqrt(delta) ) / dr^2 + p;
        I2(2,j) = (-D*dx - sgn*dy*sqrt(delta) ) / dr^2 + q;
    end
    
    % append zeros for finding true anomaly
    I1(end+1,:) = 0;
    I2(end+1,:) = 0;
    
    % find true anomaly at all intercepts
    c_vect = [p;q;0]; % vector for center of hodograph
    r1_vect = I1 - c_vect;
    r2_vect = I2 - c_vect;
    
    f1 = nan(1,size(r1_vect,2));
    for j = 1:size(r1_vect,2)
        f1(j) = atan2(norm(cross(c_vect,r1_vect(:,j))),dot(c_vect,r1_vect(:,j)));
    end
    
    f2 = nan(1,size(r2_vect,2));
    for j = 1:size(r2_vect,2)
        f2(j) = atan2(norm(cross(c_vect,r2_vect(:,j))),dot(c_vect,r2_vect(:,j)));
    end
    
    ff = [f1;f2];
    
    % get the semi-major axis. this is true regardless which value I takes
    a = mu/(I1(:,1)'*I1(:,1)) * (1 + e^2 + 2*e*cos(f1(1))) / (1 - e^2);

    % find 'real' true anomaly by propagating the two intercepts from the
    % first measurement each in two directions
    fcn = @(ee,mm) ee - e*sin(ee) - mm;
    Ea = 2 * atan(sqrt((1-e)/(1+e))*tan(f1(1)/2));
    Eb = Ea;
    Ma = Ea - e*sin(Ea);
    Mb = Ma;
    F1 = nan(2,length(dt)+1); F1(:,1) = f1(1);
    for i = 1:length(dt)
        Ma = Ma + dt(i) / sqrt(a^3/mu);
        Mb = Mb - dt(i) / sqrt(a^3/mu);
        Ea = fzero(@(ee) fcn(ee,Ma), Ma + e/2);
        Eb = fzero(@(ee) fcn(ee,Mb), Mb + e/2);
        F1(1,i+1) = 2*atan(tan(Ea/2)/sqrt((1-e)/(1+e)));
        F1(2,i+1) = 2*atan(tan(Eb/2)/sqrt((1-e)/(1+e)));
    end
    Ea = 2 * atan(sqrt((1-e)/(1+e))*tan(f2(1)/2));
    Eb = Ea;
    Ma = Ea - e*sin(Ea);
    Mb = Ma;
    F2 = nan(2,length(dt)+1); F2(:,1) = f2(1);
    for i = 1:length(dt)
        Ma = Ma + dt(i) / sqrt(a^3/mu); Ma = mod(Ma,2*pi);
        Mb = Mb - dt(i) / sqrt(a^3/mu); Mb = mod(Mb,2*pi);
        Ea = fzero(@(ee) fcn(ee,Ma), Ma + e/2);
        Eb = fzero(@(ee) fcn(ee,Mb), Mb + e/2);
        F2(1,i+1) = 2*atan(tan(Ea/2)/sqrt((1-e)/(1+e)));
        F2(2,i+1) = 2*atan(tan(Eb/2)/sqrt((1-e)/(1+e)));
    end
    
    FF = [F1;F2];
    
    % now we compare the four 'real' true anomaly paths
    % against the true anomalies calculated from the intercepts
    % and find the minimum true anomaly error
    f_err = zeros(4,length(dt));
    for i = 1:4
        for j = 1:length(dt)
            f_err(i,j) = min(abs(ff(:,j)-FF(i,j)));
        end
    end
    
    f_diff = f_err(vecnorm(f_err,2,2)==min(vecnorm(f_err,2,2)),:);
    f_diff(2:end,:) = [];
    f_diff = f_diff * 10000;
    
end
