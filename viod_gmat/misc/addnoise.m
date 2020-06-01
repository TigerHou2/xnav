function v_noisy = addnoise(v,noise,bias)
%ADDNOISE adds noise to velocity measurements provided by GMAT
%
% Author:
%   Tiger Hou
%
% Arguments:
%   v       - km/s, N-by-3 array of velocity measurements
%   noise   -  m/s, scalar value for 1 standard deviation of noise
%
% Optional Arguments:
%   bias    -  m/s, 1-by-3 array of velocity measurement bias
%
% Outputs:
%   v_noisy - km/s, N-by-3 array of perturbed velocity measurements

if ~exist('bias','var')
    bias = zeros(1,3);
end

% clean up unused velocity measurements
v_orig = v;
v = vclean(v);

noise_vec = randn(size(v));
noise_vec = noise_vec ./ vecnorm(noise_vec,2,2) * noise / 1000;

v_noisy = v + noise_vec + bias/1000;

v_noisy_resize = zeros(size(v_orig));
v_noisy_resize(1:size(v,1),1:size(v,2)) = v_noisy;
v_noisy = v_noisy_resize;

end

