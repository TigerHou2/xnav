function fVal = objFcn(vars,rangeRateData,Tvect,pulsars,mu)
%OBJFCN Objective function for range-rate OD orbit parameter guess

debug = false;

fVal = 0;

f0 = vars(1);
e = vars(2);
dM = vars(3);

if e >= 1 || e < 0
    fVal = Inf;
    return
end

period = (2*pi) / dM * ( max(Tvect(:))-min(Tvect(:)) );
a = (( period/(2*pi) )^2 * mu) ^ (1/3);
R = sqrt(mu/a/(1-e^2));

E0 = 2*atan( tan(f0/2) * sqrt((1-e)/(1+e)) );
M0 = E0 - e*sin(E0);
Mvect = ( Tvect-min(Tvect(:)) ) / period * (2*pi) + M0;
Evect = kepler(Mvect,e);
Fvect = 2*atan( tan(Evect/2) * sqrt((1+e)/(1-e)) );

numPulsars = size(Tvect,1);
pulsarCoordinates = nan(2,numPulsars);

f_pe = 0;
f_pp = 2*atan( tan(pi/2/2) * sqrt((1+e)/(1-e)) );
componentVectors = nan(2,3);

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
    fVal = fVal + norm((A*x-b)/norm(b));
    if j <= 3
        componentVectors(1,j) = R*x(1)*(e+cos(f_pe))+R*x(2)*sin(f_pe);
        componentVectors(2,j) = R*x(1)*(e+cos(f_pp))+R*x(2)*sin(f_pp);
    end
    if debug
        hold on
        angs = linspace(0,2*pi-1e-3,1000);
        plot(angs,R*x(1)*(e+cos(angs))+R*x(2)*sin(angs))
        scatter(mod(thisF,2*pi),A*x)
        hold off
    end
end

% find the orbit plane normal
v11 = pulsars(:,1) * componentVectors(1,1);
v12 = pulsars(:,2) * componentVectors(1,2);
v13 = pulsars(:,3) * componentVectors(1,3);
v21 = pulsars(:,1) * componentVectors(2,1);
v22 = pulsars(:,2) * componentVectors(2,2);
v23 = pulsars(:,3) * componentVectors(2,3);

v_pe = [v11'; v12'; v13'] \ [v11'*v11; v12'*v12; v13'*v13];
v_pp = [v21'; v22'; v23'] \ [v21'*v21; v22'*v22; v23'*v23];
% The E = 0 and E = 90 velocities are constructed in the pulsar coordinate
% system. To align the hodogrpah coordinate system to the pulsar system, we
% need to apply the convention that the periapsis velocity vector points in
% the +y direction. This arises directly from the equations that convert
% orbit parameters to position and velocity vectors.
K = cross(-v_pp, v_pe);
K = K / norm(K);
ux = -v_pp / norm(v_pp);
uy =  v_pe / norm(v_pe);

trueCoordinates = nan(2,numPulsars);
for j = 1:numPulsars
    thisPulsar = pulsars(:,j);
    thisNormal = thisPulsar'*K*K;
    thisPlanar = thisPulsar-thisNormal;
    RA  = atan2(uy'*thisPlanar,ux'*thisPlanar);
    INC = atan2(norm(thisPlanar),K'*thisNormal);
    trueCoordinates(1,j) = RA;
    trueCoordinates(2,j) = INC;
end

alignmentErrorP = trueCoordinates + pulsarCoordinates;
alignmentErrorM = trueCoordinates - pulsarCoordinates;
alignmentErrors = cat(3,alignmentErrorP,alignmentErrorM);
alignmentErrors = mod(alignmentErrors,pi/2);
alignmentErrors = cat(3,alignmentErrors,pi/2-alignmentErrors);
minAlignmentErrors = min(alignmentErrors,[],3);
AlignmentError = sum(minAlignmentErrors(:));

AngularError = 0;
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
%                 p1 = [  cos(INC_1)*cos(RA_1);...
%                         cos(INC_1)*sin(RA_1);...
%                        -sin(INC_1)];
%                 p2 = [  cos(INC_2)*cos(RA_2);...
%                         cos(INC_2)*sin(RA_2);...
%                        -sin(INC_2)];
%                 thisDiff = acos(abs(p1' * p2));
                thisDiff = acos(abs( ...
                           cos(INC_1)*cos(RA_1) * cos(INC_2)*cos(RA_2) ...
                         + cos(INC_1)*sin(RA_1) * cos(INC_2)*sin(RA_2) ...
                         + sin(INC_1) * sin(INC_2) ));
                localError = abs(thisDiff-trueDiff);
                if localErrorMin > localError
                    localErrorMin = localError;
                end
            end
%             if AngularError < localErrorMin
%                 AngularError = localErrorMin;
%             end
            AngularError = AngularError + localErrorMin;
        end
    end
end
fVal = [fVal, AlignmentError, AngularError];

end

