function expand(lm,rm,bm,um)
%EXPAND Summary of this function goes here
%   Detailed explanation goes here

if ~exist('lm','var')
    lm = 0;
end
if ~exist('rm','var')
    rm = 0;
end
if ~exist('bm','var')
    bm = 0;
end
if ~exist('um','var')
    um = 0;
end

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1) + lm;
bottom = outerpos(2) + ti(2) + bm;
ax_width = outerpos(3) - ti(1) - ti(3) - lm - rm;
ax_height = outerpos(4) - ti(2) - ti(4) - bm - um;
ax.Position = [left bottom ax_width ax_height];

end

