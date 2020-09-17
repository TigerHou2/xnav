close all
clear;clc

addpath('..\functions')
addpath('..\..\..\journal\cases')

[~,eccs,taVect,~,~] = ...
    load_orbit_cases('eccVect',linspace(0.1,0.9,17),...
                     'taVect' ,[0,135,180]);
a = 1;
mu = 1;
period = 0.1;
f0deg = 135;
f0 = deg2rad(f0deg);

% select the nth observation's position error for comparison
selObsv = 1;

period = period * 2*pi;
errVect = nan(size(eccs));

% error data for uniform TRUE anomalies

% error data for uniform MEAN anomalies
errDat = ...
{[0.00835315722751571 0.0165655659530783 0.0170780618386171;0.00695264455976074 0.0188567944992468 0.0203011089712518;0.00578243350575534 0.0212208907011933 0.0240746038627719;0.00481146299252155 0.0236025004605882 0.0284872958075576;0.00401195970850299 0.0259314554588709 0.0336466902884267;0.00335878810557408 0.0281208125877357 0.0396850320083509;0.00283149604380018 0.0300662325877522 0.0467656198199464;0.00241251433359466 0.0316454297675005 0.0550973507848616;0.00208397117458786 0.0327188881896641 0.0649558137942531;0.00183046467392142 0.0331339167507514 0.0767075857574209;0.0016402185518793 0.0327325031317418 0.0908566906527482;0.00150424459116322 0.0313661481391544 0.108130177425741;0.00141618497479626 0.0289209082567984 0.129636641164997;0.00137176419567978 0.0253570145061127 0.157183336122832;0.00136663742298531 0.0207669225572916 0.193994052906263;0.00139596954932625 0.0154478438624509 0.246646312616822;0.00145584369601623 0.00996030706316017 0.331809160765997];...
 [0.00339332423661162 0.00672990522649987 0.00693739264399444;0.00282439379246236 0.00766066854097154 0.00824656176558072;0.00234902128240103 0.00862100274014265 0.00977931693610746;0.00195458804752475 0.00958840533259112 0.0115716942623907;0.00162980442398471 0.0105343669154658 0.0136673463870168;0.00136446827905221 0.0114235619551859 0.016119955376136;0.00115026445817477 0.012213596071876 0.018995854052386;0.000980058088391366 0.0128547884620762 0.0223798578006789;0.00084659135110338 0.0132904714363779 0.0263838814138103;0.000743607396639942 0.0134586471324496 0.0311567847493022;0.000666324183724384 0.0132951778913974 0.0369032969268079;0.000611087301861676 0.012739803358355 0.0439186304630291;0.000575316641239541 0.0117463105844411 0.0526529319578112;0.000557273110191888 0.0102986075775864 0.0638399490094842;0.000555190920603286 0.00843427654203173 0.0787888523248029;0.00056710728747565 0.00627398867625942 0.100170536332225;0.000591432014279602 0.00404532590414433 0.134754327513852];...
 [0.000975299812661407 0.00193435610396258 0.00199389374556362;0.000811779516128675 0.00220187480689119 0.00237015567684784;0.000675149825199861 0.00247788916038416 0.00281067819903108;0.000561783769108566 0.00275592946960193 0.00332581688639012;0.000468435444791407 0.00302780111606707 0.00392811629124497;0.000392173779048907 0.00328335016556284 0.00463300064925917;0.000330607873659043 0.00351039057610488 0.00545953632627071;0.000281687353465874 0.00369464245468102 0.00643209611538827;0.000243326512340766 0.00381981941162268 0.00758284299185798;0.000213727002184089 0.00386810567091552 0.00895455666028608;0.000191514646070711 0.00382107304873483 0.0106060813327208;0.000175638665063574 0.00366140909182843 0.0126222530415148;0.000165357863304571 0.00337583991009489 0.0151324335616067;0.000160172080102294 0.00295974771473687 0.0183474713136107;0.000159573691244801 0.00242393742943053 0.0226436038552957;0.000162998753129845 0.00180308529541869 0.0287884076621551;0.000169990351572474 0.00116259000003698 0.0387273473056273]};

idx = find(taVect==f0);
errDat = errDat{selObsv};
if isempty(idx)
    error('The requested true anomaly value is not tabulated!')
else
    errDat = errDat(:,idx)';
end

for i = 1:length(eccs)
    e = eccs(i);
    E0 = 2 * atan(sqrt((1-e)/(1+e))*tan(f0/2));
    M0 = E0 - e*sin(E0);
    M = M0 + period;
    E = kepler(M,e);
    Mvect = linspace(M0,M,3);
    Evect = kepler(Mvect,e);
    fvect = 2 * atan(sqrt((1+e)/(1-e))*tan(Evect/2));
    fvect = mod(fvect,2*pi);
    f = fvect(end);
    df = f-f0;
    df = mod(df,2*pi);
    error = 1;
    ff = fvect(selObsv);
    % adj 1: error scales inversely with hodograph radius
    adj = sqrt((1-e)/(1+e)) + sqrt((1+e)/(1-e));
    adj = 1 / adj;
    disp(['Hodo Radius: ' num2str(adj)])
    error = error * adj;
    
%         % adj 2: error scales inversely with position magnitude
%         adj = 1 / (1-e^2) * (1+e*cos(ff));
%         disp(['Position Mag:' num2str(adj)])
%     %     error = error * adj;
%         % adj 3: error scales inversely with measurement magnitude
%         cosVal = cos(ff);
%         adj = 1 / sqrt((1-e^2) / (1+2*e*cosVal+e^2));
%         disp(['Meas Mag:    ' num2str(adj)])
%     %     error = error * adj;

    % adj 4: error scales inversely with square of measurement span
    adj = 1 / df^2;
%     adj = 1 / df^(2-df/pi);
    disp(['Meas Span:   ' num2str(adj)])
    error = error * adj;
    disp(' ')
    
    % calculate error
    errVect(i) = error;
end

xVar = eccs;
yRef = errDat;
yVar = errVect;
scaling = 1 / (max(errVect)-min(errVect)) * (max(errDat)-min(errDat));
yVar = yVar * scaling;
offset  = - min(yVar) + min(yRef);
yVar = yVar + offset;

disp(['Scaling = ' num2str(scaling)])
disp(['Offset  = ' num2str(offset)])

figure;
plot(eccs,yVar*100,'LineWidth',1.5)
hold on
plot(eccs,yRef*100,'LineWidth',1.5)
hold off
legend('Prediction','Simulation','Location','Best')
xlabel('Eccentricity')
ylabel('Position Error \%')
set(gca,'FontSize',18)
grid on