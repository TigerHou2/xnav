addpath('fcns_iod')
addpath('fcns_orb')
addpath('fcns_vis')
addpath('fcns_testing')
addpath('fcns_misc')

% astronomical unit in km
AU = 149597870;

% Sun
sun.mu = 1.327124400e11;

% Mercury
mercury.a = 0.3871;
mercury.e = 0.2056;

% Veuns
venus.a = 0.7233;
venus.e = 0.0068;

% Earth
earth.a = 1;
earth.e = 0.0167;
earth.mu = 3.986004418e5;

% Mars
mars.a = 1.5237;
mars.e = 0.0934;

% Jupiter
jupiter.a = 5.2028;
jupiter.e = 0.0483;

% Saturn
saturn.a = 9.5388;
saturn.e = 0.0560;

% Uranus
uranus.a = 19.1914;
uranus.e = 0.0461;

% Neptune
neptune.a = 30.0611;
neptune.e = 0.0097;

% Pluto
pluto.a = 39.5294;
pluto.e = 0.2482;

planets = [mercury.a, mercury.e; ...
             venus.a,   venus.e; ...
             earth.a,   earth.e; ...
              mars.a,    mars.e; ...
           jupiter.a, jupiter.e; ...
            saturn.a,  saturn.e; ...
            uranus.a,  uranus.e; ...
           neptune.a, neptune.e; ...
             pluto.a,   pluto.e];
