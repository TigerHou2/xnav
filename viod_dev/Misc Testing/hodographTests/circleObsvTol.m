% how large does noise need to be before increasing the number of
% observations actuall decreases the quality of the circle fit?
%
% we will assume angles are evenly distributed on the hodograph. This is
% because these fringe cases only happen for very short measurement times,
% so the discrepancy between mean anomaly and true anomaly is negligible.

close all
clear;clc

addpath('..\functions')

R = 1;

nvect = 0.0005:0.0005:0.01;
ovect = 3:3:100;
avect = deg2rad(1:1:18);

numSims = 1000;

dat = nan(length(nvect),1);

% iterate over all noise values
for i = 1:length(nvect)
    % iterate over all angle spans
    for j = 1:length(avect)
        errVect = nan(numSims,length(ovect));
        % iterate over all number of observations
        for k = 1:length(ovect)
            angs = linspace(0,avect(j),ovect(k))';
            obsv = R * [cos(angs) sin(angs)];
            ncube = nvect(i) * randn(ovect(k),2,numSims);
            % Monte Carlo
            for s = 1:numSims
                p = obsv + ncube(:,:,s);
                A = 2*p; A(:,3) = -1;
                B = p.^2; B = sum(B,2);
                x = A\B;
                % find radius of hodograph
                r = sqrt(x(1)^2 + x(2)^2 - x(3));
                b = x(2);
                a = x(1);
                errVect(s,k) = norm(r-R);
            end
        end
        temp = mean(errVect);
        % check if the error decreases with more observations
        if abs(max(temp)-temp(1)) < 0.05*max(temp) ...
        && abs(min(temp)-temp(length(ovect))) < 0.05 * min(temp)
            % this means the min and max are at the correct extrema
            dat(i) = rad2deg(avect(j));
            j = 1e10;
        end
    end
end
plot(nvect,dat)