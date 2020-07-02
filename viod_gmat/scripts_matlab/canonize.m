function [DUs,TUs,DoTs] = canonize(varargin)
%CANONIZE Converts to canonical units using either noise or SMA standards
%   Noise standard means the input noise corresponds to a fixed value of
%   DU/TU, whereas SMA standard means the input SMA is converted to a fixed
%   value of DU.

p = inputParser;
% reference contains 2 values:
% --- the first value is the standardization parameter in regular units
% --- the second value is the standardization parameter in canonical units
% for example, if we want to normalize a 3e-3 km/s noise to 1e-6 DU/TU, the
% --- variable will be [3e-3,1e-6].
addRequired(p,'mu')
addRequired(p,'reference')
addRequired(p,'standard')
addOptional(p,'DU',[])
addOptional(p,'TU',[])
addOptional(p,'DoT',[])
parse(p,varargin{:});

mu = p.Results.mu;
r1 = p.Results.reference(1);
r2 = p.Results.reference(2);

if strcmp(p.Results.standard,'noise')
    DU = (mu/r1^2) * r2^2;
    TU = sqrt((DU^3)/mu);
elseif strcmp(p.Results.standard,'SMA')
    DU = r1/r2;
    TU = 1/sqrt((mu/r1^3) * r2^3);
else
    warning(['No standard available for "' p.Results.standard '", '...
             'defaulting to noise standard.'])
    DU = (mu/r1^2) * r2^2;
    TU = sqrt((DU^3)/mu);
end

DUs  = p.Results.DU / DU;
TUs  = p.Results.TU / TU;
DoTs = p.Results.DoT / DU * TU;

end

