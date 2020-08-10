close all
clear;

obsVect = 3:2:150;
R = 1;
offset = [0.3 0.2];
bounds = deg2rad([0,90]);

numSims = 5000;
noise = 0.0001;

errVect = nan(numSims,length(obsVect));

for i = 1:length(obsVect)
    angs = linspace(bounds(1),bounds(2),obsVect(i))';
    dat = R * [cos(angs) sin(angs)];
    
    for s = 1:numSims
        p = dat + randn(size(dat)) * noise;
        A = 2*p; A(:,3) = -1;
        B = p.^2; B = sum(B,2);
        x = A\B;
        % find radius of hodograph
        r = sqrt(x(1)^2 + x(2)^2 - x(3));
        b = x(2);
        a = x(1);
        
        errVect(s,i) = norm(r-R);
    end
    
end

figure(1)
plot(obsVect,mean(errVect),'LineWidth',1.5)
xlabel('Number of Observations')
ylabel('Error Magnitude')
grid on

%%
f = fit(obsVect',mean(errVect)','exp2');
figure(2)
plot(f,obsVect',mean(errVect)')
f