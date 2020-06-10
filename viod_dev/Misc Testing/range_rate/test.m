close all
clear;clc

mu = 1;
a = 1;
e = 0.3;
i = 0;
omg = 0;
w = 30;

params = [a e i omg w 0];

up = [1;2;0]; % pulsar direction
up = up / norm(up);
k = [0;0;1]; % orbit normal

% calculate true values
f = deg2rad([330 0 50 200]);
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

res = 100;
rGuess = R/100;
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
        zz(i,j) = norm(fun([xx(i),yy(j),rGuess],up,v_obsv,mu,dt,k))^0.25;
    end
end

contour(X,Y,zz)
axis equal

[row,col] = find(zz==min(zz(:)));
f0 = [xx(row(1)),yy(col(1)),rGuess];

%%

f = @(x) fun(x,up,v_obsv,mu,dt,k);
options = optimoptions('fsolve','Display','iter' ...
                               ,'PlotFcn','optimplotx' ...
                               ,'MaxFunctionEvaluations',5000 ...
                               ,'MaxIterations',300);
hodo = fsolve(f,f0,options)

fun(hodo,up,v_obsv,mu,dt,k)

%% function definitions
function t_diff = fun(C,up,v,mu,dt,k)

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
            t_diff = ones(size(dt))*10000;
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
    II = cat(3,I1,I2);
    
    I = I1;
    
    % iterate through solution pairs to find ordered true anomalies
    % for now we assume measurements do not cross periapsis (f=0)

    % --- solve for true anomaly at each point
    c_vect = [p;q;0]; % vector for center of hodograph
    r1_vect = I1 - c_vect;
    r2_vect = I2 - c_vect;
    
    f1 = nan(1,size(r1_vect,2));
    for j = 1:size(r1_vect,2)
        f1(j) = atan2(norm(cross(c_vect,r1_vect(:,j))),dot(c_vect,r1_vect(:,j)));
        if dot(k,cross(c_vect,r1_vect(:,j)))<0
            f1(j) = 2*pi-f1(j);
        end
    end
    
    f2 = nan(1,size(r2_vect,2));
    for j = 1:size(r2_vect,2)
        f2(j) = atan2(norm(cross(c_vect,r2_vect(:,j))),dot(c_vect,r2_vect(:,j)));
        if dot(k,cross(c_vect,r2_vect(:,j)))<0
            f2(j) = 2*pi-f2(j);
        end
    end
    
    ff = [f1;f2];
    
    % --- find ordered true anomaly
    indices = dec2bin(0:1:2^length(v)-1)-'0' + 1;
    f_sorted = [];
    idx_sorted = [];
    for i = 1:size(indices,1)
        idx = sub2ind(size(ff),indices(i,:),1:length(v));
        if length(longestMono([ff(idx),ff(idx)])) == length(dt)+1
            f_sorted(end+1,:) = ff(idx);
            idx_sorted(end+1,:) = indices(i,:);
        end
%         if issorted(ff(idx),'monotonic')
%             f_sorted(end+1,:) = ff(idx);
%             idx_sorted(end+1,:) = indices(i,:);
%         end
    end
    
    % --- get rid of rows with duplicate values
    f_tol = 1e-3;
    idx_sorted(sum(abs(diff(sort(f_sorted,2),1,2))>f_tol,2)+1~=length(v),:) = [];
    f_sorted(sum(abs(diff(sort(f_sorted,2),1,2))>f_tol,2)+1~=length(v),:) = [];
    
    % --- if no sorted true anomalies exist, pick a random one
    if isempty(f_sorted)
        f = f1;
        I = I1;
        
        a = mu/(I(:,1)'*I(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);
    
        E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
        M = E - e*sin(E);
        t = sqrt(a^3/mu)*M;
        dt_guess = t(2:end) - t(1:end-1);
        dt_guess = mod(dt_guess,sqrt(a^3/mu));

        t_diff = dt_guess - dt;
        
        return
    end
    
    % --- otherwise, iterate through all sorted sets and find minimum time
    t_diff_vect = [];
    for j = 1:size(f_sorted,1)
        f = f_sorted(j,:);
        index = idx_sorted(j,:);
        for i = 1:length(v)
            I(:,i) = II(:,i,index(i));
        end
        
        a = mu/(I(:,1)'*I(:,1)) * (1 + e^2 + 2*e*cos(f(1))) / (1 - e^2);
    
        E = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
        M = E - e*sin(E);
        t = sqrt(a^3/mu)*M;
        dt_guess = t(2:end) - t(1:end-1);

        t_diff_vect(end+1,:) = dt_guess - dt;
    end

    t_diff = t_diff_vect(vecnorm(t_diff_vect,2,2)==min(vecnorm(t_diff_vect,2,2)),:);
    
end

function V = longestMono(S)
V = [];
for k = 1:numel(S)
    recfun(S(k),S(k+1:end))
end
% Recursive function:
    function recfun(Z,S)
    if numel(Z)>numel(V)
        V = Z;
    end
    for k = 1:numel(S)
        if Z(end)<S(k)
            recfun([Z,S(k)],S(k+1:end))
        end
    end
    end
end
