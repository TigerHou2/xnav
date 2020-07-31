close all
clear;clc

noise = 0.00005;

R = [1 5 10 20];
ang = [-5 0 5];
ang = deg2rad(ang);

C = [-50;50];

errVect = [];
cErrVect = [];

for i = 1:length(R)
    RR = R(i);
    error = 0;
    cError = 0;
    numSims = 10000;
    for n = 1:numSims
        dat = RR * [cos(ang); sin(ang)] + C;
        dat = dat + randn(size(dat)) * noise;
        % alternative circle fitting method
        A = 2*dat'; A(:,3) = -1;
        B = dat'.^2; B = sum(B,2);
        x = A\B;
        % find radius of hodograph
        r = sqrt(x(1)^2 + x(2)^2 - x(3));
        a = x(1);
        b = x(2);
        
        error = error + abs(RR-r);
        cError = cError + norm([a;b]-C);
    end
    error = error / numSims;
    cError = cError / numSims;
    errVect(end+1) = error;
    cErrVect(end+1) = cError;
end
errVect
cErrVect