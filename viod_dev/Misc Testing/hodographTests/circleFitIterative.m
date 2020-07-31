%%
% tests circle fitting by iteratively centering the circle

close all
clear;clc

R = 1;
f0 = deg2rad(30);

duration = 0.1;
numObsv = 6;
numSims = 3000;

n = 2.917e-3;
% n = 1e-4;

e = 0.4;
C = [1 1] / norm([1 1]) * R * e;

E = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
M = E - e*sin(E);
% determine measurement locations
Mvect = M + linspace(0,duration,numObsv) * 2*pi;
Evect = kepler(Mvect,e);
fVect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));

Re = 0;
Rc = 0;

diff = 1;
offset = [0 0];
itr = 1;

noise = randn(length(fVect),2) * n;

while diff > 1e-3
    
    P = R * [sin(fVect') cos(fVect')] + C - offset;
    P = P + noise;

    % construct linear system [1] eqn.9
    A = 2*P;
    A(:,3) = -1;
    B = P.^2;
    B = sum(B,2);
    x = A\B;

    Rg = sqrt(x(1)^2 + x(2)^2 - x(3)); % [1] eqn.10
    Cg = [x(1), x(2)]; % [1] eqn.11

    Re = abs(R-Rg);
    Rc = norm(C-Cg);
    
    offset = offset + Cg;
    diff = norm(Cg);
    
    disp(['itr = ' num2str(itr) ', Re = ' num2str(Re)])
    disp(Cg)
    disp(' ')
    itr = itr + 1;
    
end



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