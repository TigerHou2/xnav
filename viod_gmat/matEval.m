function out = matEval(command,varargin)
%MATEVAL Gives GMAT access to the MATLAB suite of operators and functions.
%
% Author:
%   Tiger Hou
%

for index_ = 1:nargin-1
    assignin('caller',inputname(index_+1),varargin{index_});
end

out = evalin('caller',command);

end

