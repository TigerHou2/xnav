%% nbody_test_cases
% Sources: 
%   - 13 stable 3-body cases (cases 1-13)
%       https://arxiv.org/pdf/1303.0181v1.pdf
%   - figure eight (case 18)
%       https://www.jstor.org/stable/2661357?seq=1
%       https://en.wikipedia.org/wiki/Three-body_problem#cite_note-11
%   - Apollo 11 Free Return (case 19)
%       https://space.stackexchange.com/questions/20151/lunar-free-return-trajectory-simulation
%
% other interesting cases to explore:
%   - https://www.ams.org/notices/200105/fea-montgomery.pdf
%   - http://www.maia.ub.edu/dsg/2001/index.html (Carles Simo, LOTS of n-body examples)
%
%
% from none other than the Flat Earth Society:
%   - https://wiki.tfes.org/Three_Body_Problem

%% Define Defaults & Initialize

G = 1;
r0 = zeros(3,2);
v0 = zeros(3,2);
m = [1,1,1];
tk = 0; % starting time
T = 500; % final time
dt = 1e-4; % time step

%% Define Cases

switch case_select

    case 1
        % ========= Butterfly I =========
        v0(1,1) = 0.30689;
        v0(2,1) = 0.12551;
        T = 6.2356;
        dt = T / 200000;

    case 2
        % ========= Butterfly II =========
        v0(1,1) = 0.39295;
        v0(2,1) = 0.09758;
        T = 7.0039;
        dt = T / 500000;

    case 3
        % ========= Bumblebee =========
        v0(1,1) = 0.18428;
        v0(2,1) = 0.58719;
        T = 63.5345;
        dt = T / 1000000;
    
    case 4
        % ========= Moth I =========
        v0(1,1) = 0.46444;
        v0(2,1) = 0.39606;
        T = 14.8939;
        dt = T / 100000;
        
    case 5
        % ========= Moth II =========
        v0(1,1) = 0.43917;
        v0(2,1) = 0.45297;
        T = 28.6703;
        dt = T / 100000;
        
    case 6
        % ========= Butterfly III =========
        v0(1,1) = 0.40592;
        v0(2,1) = 0.23016;
        T = 13.8658;
        dt = T / 1000000;
        
    case 7
        % ========= Moth III =========
        v0(1,1) = 0.38344;
        v0(2,1) = 0.37736;
        T = 25.8406;
        dt = T / 1000000;
        
    case 8
        % ========= Goggles =========
        v0(1,1) = 0.08330;
        v0(2,1) = 0.12789;
        T = 10.4668;
        dt = T / 100000;
        
    case 9
        % ========= Butterfly IV =========
        v0(1,1) = 0.350112;
        v0(2,1) = 0.07934;
        T = 79.4759;
        dt = T / 1000000;
        
    case 10
        % ========= Dragonfly =========
        v0(1,1) = 0.08058;
        v0(2,1) = 0.58884;
        T = 21.2710;
        dt = T / 1000000;
        
    case 11
        % ========= Yarn =========
        v0(1,1) = 0.55906;
        v0(2,1) = 0.34919;
        T = 55.5018;
        dt = T / 2000000;
        
    case 12
        % ========= Yin-Yang Ia =========
        v0(1,1) = 0.51394;
        v0(2,1) = 0.30474;
        T = 17.3284;
        dt = T / 500000;
        
    case 13
        % ========= Yin-Yang IIa =========
        v0(1,1) = 0.41682;
        v0(2,1) = 0.33033;
        T = 55.7898;
        dt = T / 1000000;
        
    case 14
        % ======== Two-Body ========
        r0 = [  0,-50, 0; ...
            0, 50, 0]';
        v0 = [-2, 0, 0; ...
            2, 0, 0]';
        m = [1,1];
        T = 1500;
        dt = 1e-1;
        G = 1000;
        
    case 15
        % ======== P-Type Binary ========
        r0 = [  0,-50, 0; ...
                0, 50, 0; ...
                0, 375, 0]';
        v0 = [ -2, 0, 0; ...
                2, 0, 0; ...
                2, 0, 0]';
        m = [1,1,0.01];
        T = 30000;
        dt = 1e-1;
        G = 1000;
        
    case 16
        % ======== S-Type Binary ========
        r0 = [  0,-50, 0; ...
                0, 50, 0; ...
                0, 47, 0]';
        v0 = [ -2, 0, 0; ...
                2, 0, 0; ...
               20, 0, 0]';
        m = [1,1,0.01];
        T = 500;
        dt = 1e-3;
        G = 1000;
        
    case 17
        % ======== Orbit Resonance ========
        r0 = [  0, 0, 0; ...
                0, 150, 0; ...
                0,-114, 0]';
        v0 = [ 0, 0, 0; ...
             -20, 0, 0; ...
            17.5, 0, 0]';
        m = [50,1,1];
        T = 700;
        dt = 4e-2;
        G = 1000;
        
    case 18
        % ======== Figure Eight ========
        r0(1,1) = -0.97000436;
        r0(1,2) =  0;
        r0(1,3) =  0.97000436;
        r0(2,1) =  0.24308753;
        r0(2,2) =  0;
        r0(2,3) = -0.24308753;
        v0(1,1) =  0.4662036850;
        v0(1,2) = -0.93240737;
        v0(1,3) =  0.4662036850;
        v0(2,1) =  0.4323657300;
        v0(2,2) = -0.86473146;
        v0(2,3) =  0.4323657300;
        
        m = [1,1,1];
        T = 2*pi;
        dt = T / 100000;
        G = 1;
        
    case 19
        % ======== Apollo 11 Free Return ========
        r0 = [0.4240363252016235, -0.9248449798862485, -1.232690294681233e-4;...
              0.4220528422463315, -0.9230209264977778,  1.632323615688905e-5;...
              0.4240447232851519, -0.9247715402118077, -1.129301018611092e-4]';
        v0 = [5.622675894714279, 2.5745894556521574, 3.8057228235271535e-4;...
              5.486589374929882, 2.420601498441581, -0.014677846271227611;...
              4.395253850175561, 3.8323649107803948, 0.15792573886687206]';
        m = [0.000003003, 3.69396868e-8, 5.9747e-22];
        G = 39.5;
        T = 0.01655;
        dt = T / 100000;
        lgd = {'Earth','Moon','Apollo 11'};
        
    case 20
        nbody_test_cases_extended;
    
    otherwise
        disp('Option invalid, defaulting to Butterfly I.')
        % ========= Butterfly I =========
        v0(1,1) = 0.30689;
        v0(2,1) = 0.12551;
        T = 6.2356;
        
end

%% Finalize Initial Conditions

% fill in remaining initial conditions for first 13 cases
if case_select < 14
    r0(1,1) = -1;
    r0(1,2) =  1;
    r0(1,3) =  0;
    v0(1,2) =  1 * v0(1,1);
    v0(1,3) = -2 * v0(1,1);
    v0(2,2) =  1 * v0(2,1);
    v0(2,3) = -2 * v0(2,1);
end

% if no legend is provided, default to particle naming scheme
if ~logical(exist('lgd','var'))
    lgd{size(r0,2)} = '';
    for i = 1:size(r0,2)
        lgd{i} = strcat('Particle ', num2str(i));
    end
end