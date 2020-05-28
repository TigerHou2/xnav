function scatter3_gmat(name)
%SCATTER3_GMAT scatter3() function designed to be called from GMAT
%
% Author:
%   Tiger Hou

temp = evalin('base',name);
temp = vclean(temp);

figure
scatter3(temp(:,1),temp(:,2),temp(:,3))

end

