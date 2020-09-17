function error = viodSim(mu,a,e,i,o,w,f,dm,noise,numObsv,numSims)
%VIODSIM Summary of this function goes here
%   Detailed explanation goes here

orbitParams = [a,w,i,o,w,f];
v = nan(numObsv,3);
ncube = randn(numObsv,3,numSims);
ncube = ncube ./ vecnorm(ncube,2,2) * noise;
errDat = nan(1,numSims);
selObsv = 1;

% find measurement positions by true anomaly
E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f/2));
M0 = E0 - e*sin(E0);
M = M0 + dm;
Mvect = linspace(M0,M,numObsv);
Evect = kepler(Mvect,e);
fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
fvect = mod(fvect,2*pi);

% ground truth data
for j = 1:numObsv
    orbitParams(6) = fvect(j);
    [~,v(j,:)] = Get_Orb_Vects(orbitParams,mu);
end

% get position reference
orbitParams(6) = fvect(selObsv);
rRef = Get_Orb_Vects(orbitParams,mu);

% Monte Carlo
for s = 1:numSims
    nvect = ncube(:,:,s);
    r = hodo(v+nvect,mu);
    r = r(selObsv,:)';
    errDat(s) = norm(r-rRef) / norm(rRef);
end

error = mean(errDat);

end

