function append_data(name,data,dim)
%APPEND_DATA Appends data to a workspace variable specified in GMAT
%
% Author:
%   Tiger Hou
%
% Description:
%   This function is designed to be called from GMAT. 
%   The user specifies a workspace variable using a string, a set of data
%   to append to that workspace variable, and optionally, the dimension
%   along which the data should be appended.

if ~exist(name,'var')
    assignin('base',name,[]);
end

if ~exist('dim','var')
    dim = 1;
end

temp = evalin('base',name);
temp = cat(dim,temp,data);
assignin('base',name,temp);

end

