function v = vclean(v)
%VCLEAN Clean up extra entries in velocity measurements
%
% Author:
%   Tiger Hou
%
% Description:
%   This function is necessary because GMAT doesn't allow dynamic variable
%   sizing. It also does not allow us to specify the variable size using a
%   variable. 
%   The solution is to preallocation 100 velocity measurements in GMAT.
%   Then, we use this function to clean up any extra values we do not need.

v(sum(v ~= 0,2)==0,:) = [];

end

