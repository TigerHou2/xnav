function fVal = objFcn(vars,rangeRateData,Tvect,pulsars,mu)
%OBJFCN Objective function for range-rate OD orbit parameter guess

fVal = 0;

f0 = vars(1);
e = vars(2);
dM = vars(3);

period = (2*pi) / dM * ( max(Tvect(:))-min(Tvect(:)) );
a = (( period/(2*pi) )^2 * mu) ^ (1/3);
R = sqrt(mu/a/(1-e^2));

E0 = 2*atan( tan(f0/2) * sqrt((1-e)/(1+e)) );
M0 = E0 - e*sin(E0);
Mvect = ( Tvect-min(Tvect(:)) ) / period * (2*pi) + M0;
Evect = kepler(Mvect,e);
Fvect = 2*atan( tan(Evect/2) * sqrt((1+e)/(1-e)) );

numPulsars = size(Tvect,1);
guessPulsars = nan(size(pulsars));
pulsarCoordinates = nan(2,numPulsars);

for j = 1:numPulsars
    thisF = Fvect(j,:);
    A = R*[e+cos(thisF)', sin(thisF)'];
    b = rangeRateData(j,:)';
    x = A\b;
    x(1) = max(-1,min(1,x(1)));
    x(2) = max(-1,min(1,x(2)));
    RA = atan(x(1)/x(2));
    if abs(cos(RA)) < 0.05
        INC = acos(x(1)/sin(RA));
    else
        INC = acos(x(2)/cos(RA));
    end
    pulsarCoordinates(:,j) = real([RA;INC]);
    fVal = fVal + norm(A*x-b);
end

maxAngularError = 0;
pm = [  1,  1,  1,  1; ... 16 combinations, only 4 unique cases
        1,  1,  1, -1; ... 
        1, -1,  1,  1; ... 
        1, -1,  1, -1; ... 
      ]';
for j = 1:numPulsars
    for k = 1:numPulsars
        if j>k
            localErrorMin = 1;
            trueDiff = acos(abs(pulsars(:,j)' * pulsars(:,k)));
            for i = 1:4
                RA_1  = pulsarCoordinates(1,j) * pm(1,i);
                RA_2  = pulsarCoordinates(1,k) * pm(2,i);
                INC_1 = pulsarCoordinates(2,j) * pm(3,i);
                INC_2 = pulsarCoordinates(2,k) * pm(4,i);
                p1 = [  cos(INC_1)*cos(RA_1);...
                        cos(INC_1)*sin(RA_1);...
                       -sin(INC_1)];
                p2 = [  cos(INC_2)*cos(RA_2);...
                        cos(INC_2)*sin(RA_2);...
                       -sin(INC_2)];
                thisDiff = acos(abs(p1' * p2));
                localError = abs(thisDiff-trueDiff);
                if localErrorMin > localError
                    localErrorMin = localError;
                end
            end
            if maxAngularError < localErrorMin
                maxAngularError = localErrorMin;
            end
        end
    end
end
fVal = fVal * maxAngularError;

end

