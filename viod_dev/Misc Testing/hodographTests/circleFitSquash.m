% compares circle fitting algorithm accuracy w.r.t separation of samples

close all
clear;

addpath('functions')

% circle definition
R = 1;
e = 0.1;
C = e * [1;1] / norm([1;1]) * R;

% define ranges
range = logspace(log10(8),log10(180),150);
% range = linspace(8,180,150);
errVect = nan(size(range));

sample_count = 3;
f0 = 90;

while true

% iterate through ranges
for j = 1:length(range)

% take samples
sample_range = [0 range(j)] + f0;

    sample_range = deg2rad(sample_range);
    E_range = 2 * atan(sqrt((1-e)/(1+e))*tan(sample_range/2));
    M_range = E_range - e*sin(E_range);
    M_vect = linspace(M_range(1),M_range(2),sample_count);
    E_vect = kepler(M_vect,e);
    f_vect = 2 * atan(sqrt((1+e)/(1-e))*tan(E_vect/2));
    
    f_vect = linspace(sample_range(1),sample_range(2),sample_count);

points = [cos(f_vect);sin(f_vect)] + C;

err = 0;

% Monte Carlo simulation
numSims = 750;
for i = 1:numSims
    
    % apply noise
    noise = 0.0005;
    p = points + randn(size(points))*noise;

%     circle fitting
%     [a,b,r] = hyperfit(p);
    
        % alternative circle fitting method
        A = 2*p'; A(:,3) = -1;
        B = p'.^2; B = sum(B,2);
        x = A\B;
        % find radius of hodograph
        r = sqrt(x(1)^2 + x(2)^2 - x(3));
        b = x(2);
        a = x(1);
    
    % calculate error
    err = err + abs(R-r);
%     err = err + norm(C-[a;b]);
    
%     if r > 100*R
%         warning(['Radius estimate is very bad! (r = ' num2str(r) ')'])
%     end
    
end

errVect(j) = err / numSims;

end

% does the inverse-square model fit?
x = 1./(range'.^2);
A = [x, ones(size(x))];
y = errVect';
B = y;

ymean = mean(y);
sstot = sum((y-ymean).^2);
mb = A\B;
ssres = sum((A*mb-B).^2);
R2 = 1-ssres/sstot;

% plot
h = plot(x,y);
xlabel('$1/df^2$');
ylabel('error');
title(['f0 = ' num2str(f0) ',   R = ' num2str(R2,4)])
drawnow;

f0 = mod(f0+3,360);
break
if ~isvalid(h)
    break
end

end

% display line properties



