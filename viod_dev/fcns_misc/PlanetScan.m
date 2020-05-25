function planet = PlanetScan(r1,r2,planets)
%PLANETSCAN Returns the names of planets having r1 and r2 semi-major axes
%   planets is a 9x2 matrix where the first column contains semi-major axis
%   information for the eight planets and Pluto, and the second column
%   contains corresponding eccentricity info.

planet.start = '';
planet.end = '';

p1 = find(planets==r1, 1);
p2 = find(planets==r2, 1);

if isempty(p1) || isempty(p2)
    planet = [];
    return
end

names = {'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', ...
         'Saturn', 'Uranus', 'Neptune', 'Pluto'};

planet.start = names{p1};
planet.end   = names{p2};

end

