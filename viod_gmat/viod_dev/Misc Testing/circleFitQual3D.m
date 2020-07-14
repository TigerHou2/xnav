%%
% tests circle fitting quality dependence on separation of sample points
% this is done by varying eccentricity

close all
clear;clc

R = 1;

numObsv = 3;
numSims = 1000;

n = 2.917e-3;
% n = 1e-4;

f0vect  = deg2rad(linspace(0,180,20));
durVect = linspace(0.1,0.9,20);

eccs = linspace(0,0.95,20);
ReVect = nan(size(eccs));
RcVect = nan(size(eccs));

f0Dat  = nan(length(f0vect)*length(durVect)*length(eccs),1);
durDat = nan(length(f0vect)*length(durVect)*length(eccs),1);
eccDat = nan(length(f0vect)*length(durVect)*length(eccs),1);
errDat = nan(length(f0vect)*length(durVect)*length(eccs),1);
stdDat = nan(length(f0vect)*length(durVect)*length(eccs),1);

combVec = combvec(f0vect,durVect,eccs);
combIdx = combvec(1:20,1:20,1:20);

parfor q = 1:size(combVec,2)
    i = combIdx(1,q);
    j = combIdx(2,q);
    k = combIdx(3,q);
    f0 = f0vect(i);
    duration = durVect(j);
    e = eccs(k);
    
C = [1 1] / norm([1 1]) * R * e;

E = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M = E - e*sin(E);
% determine measurement locations
Mvect = M + linspace(0,duration,numObsv) * 2*pi;
Evect = kepler(Mvect,e);
fVect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));

Re = nan(numSims,1);
Rc = nan(numSims,1);

for m = 1:numSims

    P = R * [sin(fVect') cos(fVect')] + C;
    P = P + randn(size(P)) * n;

    % construct linear system [1] eqn.9
    A = 2*P;
    A(:,3) = -1;
    B = P.^2;
    B = sum(B,2);
    x = A\B;

    Rg = sqrt(x(1)^2 + x(2)^2 - x(3)); % [1] eqn.10
    Cg = [x(1), x(2)]; % [1] eqn.11
    
    Re(m) = abs(R-Rg);
    Rc(m) = norm(C-Cg);
    
end

f0Dat(q)  = f0;
durDat(q) = duration;
eccDat(q) = e;
errDat(q) = mean(Re);
stdDat(q) = std(Re);

end

%% plot results

figure;
colormap hsv
% size = (log(errDat) - min(log(errDat)) + 1).^1;
% size = size / max(size) * 30;
pos = 0;
tol = 0.2;
lb = (max(log(errDat))-min(log(errDat)))*(pos-tol) + min(log(errDat));
ub = (max(log(errDat))-min(log(errDat)))*(pos+tol) + min(log(errDat));
ff = rad2deg(f0Dat);
dd = durDat;
ee = eccDat;
rr = log(errDat);
ff(rr>ub|rr<lb) = [];
dd(rr>ub|rr<lb) = [];
ee(rr>ub|rr<lb) = [];
rr(rr>ub|rr<lb) = [];
ptSize = 12;
fig = scatter3(ff,dd,ee,ptSize,rr,'filled');
fig.MarkerFaceAlpha = 1;
fig.MarkerEdgeAlpha = 1;
xlabel('f0')
ylabel('duration')
zlabel('eccentricity')
colorbar
pbaspect([1 1 1])

viewDir = [1 1 1]';
ang = 1;
rot = [cosd(ang) -sind(ang) 0;...
       sind(ang)  cosd(ang) 0;...
       0          0         1];

while true
    viewDir = rot * viewDir;
    view(viewDir)
    drawnow
    if ~isvalid(fig)
        break
    end
    pause(1/120)
end

clf
colormap hsv
fig = scatter3(ff,dd,ee,ptSize,rr,'filled');
fig.MarkerFaceAlpha = 1;
fig.MarkerEdgeAlpha = 1;
xlabel('f0')
ylabel('duration')
zlabel('eccentricity')
colorbar
pbaspect([1 1 1])


%% function definitions
function Evect = kepler(Mvect,e)
%KEPLER solves Kepler's equation using the Conway-Laguerre algorithm.

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