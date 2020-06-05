function [R1,V1] = predict_ecc(V,mu,noise,R,P,E,ecc_idx)
%PREDICT_ECC Plots error statistics for VIOD w.r.t eccentricity using a
%GMAT orbit propagator.
%
% Author:
%   Tiger Hou
%

% declare persistent variables
persistent rms per ecc ecc_vect

R1 = nan;
V1 = nan;

%% OD mode
if nargin == 3
    
    % clean up unused velocity measurement slots
    V = vclean(V);
    
    vn = addnoise(V,noise);
    rn = viod(vn,mu);
    R1 = rn(1,:);
    V1 = vn(1,:);
    
    return;
    
%% data storage mode
elseif nargin == 7
    
    rms(end+1,:) = R;
    per(end+1,:) = P;
    ecc(end+1,:) = E;
    if ~any(ismember(ecc_vect,ecc_idx))
        ecc_vect(1,end+1) = ecc_idx;
    end
    
    return;
    
%% data plotting mode
elseif nargin == 0
    
    ecc_mat = repmat(ecc_vect,numel(ecc)/length(ecc_vect),1);
    
    rms = vecnorm(rms,2,2);
    rms = reshape(rms,numel(rms)/length(ecc_vect),length(ecc_vect));
    per = reshape(per,numel(per)/length(ecc_vect),length(ecc_vect));
    ecc = reshape(ecc,numel(ecc)/length(ecc_vect),length(ecc_vect));

    linewidth = 1;
    scFormat = 'r.';
    errFormat = '-k+';

    figure;
    latexify(45,18)

    subplot(1,3,1)
    hold on
    scatter(ecc_mat(:),rms(:),scFormat)
    errorbar(ecc_vect,mean(rms),std(rms),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Eccentricity, nd')
    ylabel('RMS Position Error, $km$')

    subplot(1,3,2)
    hold on
    scatter(ecc_mat(:),per(:),scFormat)
    errorbar(ecc_vect,mean(per),std(per),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Eccentricity, nd')
    ylabel('Period Error, seconds')

    subplot(1,3,3)
    hold on
    scatter(ecc_mat(:),ecc(:),scFormat)
    errorbar(ecc_vect,mean(ecc),std(ecc),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Eccentricity, nd')
    ylabel('Eccentricity Error, nd')
    
    % clean up persistent variables
    clear rms per ecc ecc_vect

end

end
