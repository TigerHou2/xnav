function [r,v,initGuess] = rrod(rangeRateData,Tvect,pulsars,mu)
%RROD Performs range-rate orbit determination.
%   Author: Tiger Hou

fun = @(x) obj_fcn(x,rangeRateData,Tvect,pulsars,mu);
options = optimset( 'Display','none', ...
                    'TolFun',1e-16, 'TolX', 1e-10);
                
initGuess = [pi, 0.5, pi];
radius = [pi+0.1, 0.5, pi];

guesses = nan(3,3);
for ii = 1:3
resolution = [35,35,35];
initGuess(ii) = initGuess(ii) + radius(ii)/resolution(ii);
lb = initGuess - radius;
ub = initGuess + radius;
range_f0 = linspace(lb(1),ub(1),resolution(1)+1);
range_e  = linspace(lb(2),ub(2),resolution(2)+1);
range_dM = linspace(lb(3),ub(3),resolution(3)+1);

temp = local_search(fun,range_f0,range_e,range_dM,resolution);
temp = fminsearch(fun,temp,options);
temp(1) = mod(temp(1),2*pi);

guesses(ii,:) = temp;
end
[~,idx] = min([fun(guesses(1,:)),fun(guesses(2,:)),fun(guesses(3,:))]);
initGuess = guesses(idx,:);

initGuess(1) = mod(initGuess(1),2*pi);
[rot,params] = get_rot(initGuess,rangeRateData,Tvect,pulsars,mu);
[r,v] = Get_Orb_Vects(params,mu);
r = rot'*r;
v = rot'*v;

end

%%
function [initGuess,fVal] = ...
    local_search(fun,f0Vect,eVect,dMVect,res)
%LOCAL_SEARCH searches a given space for fmin at a fixed resolution.

fVal = Inf;
initGuess = [pi,0.5,pi];
for i = 1:res(1)
for j = 1:res(2)
for k = 1:res(3)
    in = [f0Vect(i),eVect(j),dMVect(k)];
    fValThis = fun(in);
    if fValThis < fVal
        fVal = fValThis;
        initGuess = in;
    end
end
end
end

end

%%
function fVal = obj_fcn(guess,rangeRateData,Tvect,pulsars,mu)
%OBJFCN Objective function for range-rate OD orbit parameter guess

fVal = 0;
f0 = guess(1);
e = guess(2);
dM = guess(3);

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
    thisF = Fvect(j,:)';
    A = R*[e+cos(thisF), sin(thisF)];
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
end

% find the orbit plane normal
v11 = pulsars(:,1)' * componentVectors(1,1);
v12 = pulsars(:,2)' * componentVectors(1,2);
v13 = pulsars(:,3)' * componentVectors(1,3);
v21 = pulsars(:,1)' * componentVectors(2,1);
v22 = pulsars(:,2)' * componentVectors(2,2);
v23 = pulsars(:,3)' * componentVectors(2,3);

v_pe = [v11; v12; v13] \ [v11*v11'; v12*v12'; v13*v13'];
v_pp = [v21; v22; v23] \ [v21*v21'; v22*v22'; v23*v23'];
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
            trueDiff = pulsars(:,j)' * pulsars(:,k);
            for i = 1:4
                RA_1  = pulsarCoordinates(1,j) * pm(1,i);
                RA_2  = pulsarCoordinates(1,k) * pm(2,i);
                INC_1 = pulsarCoordinates(2,j) * pm(3,i);
                INC_2 = pulsarCoordinates(2,k) * pm(4,i);
                thisDiff = cos(INC_1)*cos(RA_1) * cos(INC_2)*cos(RA_2) ...
                         + cos(INC_1)*sin(RA_1) * cos(INC_2)*sin(RA_2) ...
                         + sin(INC_1) * sin(INC_2);
                localError = abs(thisDiff-trueDiff);
                if localErrorMin > localError
                    localErrorMin = localError;
                end
            end
            AngularError = AngularError + localErrorMin;
        end
    end
end
fVal = fVal + AlignmentError + AngularError;

end

%%
function [r,v] = Get_Orb_Vects(params,mu)
%GET_ORB_VECTS Finds position and velocity vectors from orbital parameters.

a = params(1);
e = params(2);
i = params(3);
omg = params(4);
w = params(5);
f = params(6);

rr = a*(1-e^2)/(1+e*cos(f));
h = sqrt(mu*a*(1-e^2));

r = rr * [ cos(omg)*cos(w+f) - sin(omg)*sin(w+f)*cos(i) ;...
           sin(omg)*cos(w+f) + cos(omg)*sin(w+f)*cos(i) ;...
           sin(w+f)*sin(i) ];
v = -mu/h * ...
    [ cos(omg)*(sin(w+f)+e*sin(w)) + sin(omg)*(cos(w+f)+e*cos(w))*cos(i) ;...
      sin(omg)*(sin(w+f)+e*sin(w)) - cos(omg)*(cos(w+f)+e*cos(w))*cos(i) ;...
    -(cos(w+f)+e*cos(w))*sin(i) ];

end

%%
function [rot,params] = get_rot(guess,rangeRateData,Tvect,pulsars,mu)
%GET_ROT returns the orbit orientation relative to the pulsar frame

f0 = guess(1);
e = guess(2);
dM = guess(3);

period = (2*pi) / dM * ( max(Tvect(:))-min(Tvect(:)) );
a = (( period/(2*pi) )^2 * mu) ^ (1/3);
R = sqrt(mu/a/(1-e^2));

E0 = 2*atan( tan(f0/2) * sqrt((1-e)/(1+e)) );
M0 = E0 - e*sin(E0);
Mvect = ( Tvect-min(Tvect(:)) ) / period * (2*pi) + M0;
Evect = kepler(Mvect,e);
Fvect = 2*atan( tan(Evect/2) * sqrt((1+e)/(1-e)) );

numPulsars = size(Tvect,1);

f_pe = 0;
f_pp = 2*atan( tan(pi/2/2) * sqrt((1+e)/(1-e)) );
componentVectors = nan(2,3);

for j = 1:numPulsars
    thisF = Fvect(j,:)';
    A = R*[e+cos(thisF), sin(thisF)];
    b = rangeRateData(j,:)';
    x = A\b;
    x(1) = max(-1,min(1,x(1)));
    x(2) = max(-1,min(1,x(2)));
    if j <= 3
        componentVectors(1,j) = R*x(1)*(e+cos(f_pe))+R*x(2)*sin(f_pe);
        componentVectors(2,j) = R*x(1)*(e+cos(f_pp))+R*x(2)*sin(f_pp);
    end
end

% find the orbit plane normal
v11 = pulsars(:,1)' * componentVectors(1,1);
v12 = pulsars(:,2)' * componentVectors(1,2);
v13 = pulsars(:,3)' * componentVectors(1,3);
v21 = pulsars(:,1)' * componentVectors(2,1);
v22 = pulsars(:,2)' * componentVectors(2,2);
v23 = pulsars(:,3)' * componentVectors(2,3);

v_pe = [v11; v12; v13] \ [v11*v11'; v12*v12'; v13*v13'];
v_pp = [v21; v22; v23] \ [v21*v21'; v22*v22'; v23*v23'];
% The E = 0 and E = 90 velocities are constructed in the pulsar coordinate
% system. To align the hodogrpah coordinate system to the pulsar system, we
% need to apply the convention that the periapsis velocity vector points in
% the +y direction. This arises directly from the equations that convert
% orbit parameters to position and velocity vectors.
K = cross(-v_pp, v_pe);
K = K / norm(K);
ux = -v_pp / norm(v_pp);
uy =  v_pe / norm(v_pe);

rot = [ux, uy, K]';
params = [a,e,0,0,0,f0];

end
