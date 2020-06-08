function var_ecc(V,mu,noise,num_monte,pos_ref,period_ref,ecc_ref)
%VAR_ECC Plots error statistics for VIOD w.r.t eccentricity.
%
% Author:
%   Tiger Hou
%

% declare persistent variables
persistent rms period ecc ecc_mat

%% data collection mode
if nargin == 7
    
    rms_temp = nan(num_monte,1);
    period_temp = nan(num_monte,1);
    ecc_temp = nan(num_monte,1);
    
    % clean up unused velocity measurement slots
    V = vclean(V);
    
    for i = 1:num_monte

        vn = addnoise(V,noise);
        rn = viod(vn,mu);
        rr = rn(1,:);
        vv = vn(1,:);

        rms_temp(i)    = norm(rr-pos_ref);
        sma_est        = norm(rr)/(2 - norm(rr)*dot(vv,vv)/mu);
        period_temp(i) = 2*pi*sqrt(sma_est^3/mu) - period_ref;
        ecc_temp(i)    = norm(((dot(vv,vv)-mu/norm(rr))*rr-dot(rr,vv)*vv)/mu)...
                       - ecc_ref;

    end
    
    rms = [rms rms_temp];
    period = [period period_temp];
    ecc = [ecc ecc_temp];
    ecc_mat = [ecc_mat repmat(ecc_ref,num_monte,1)];
    
    return;
    
%% data plotting mode
else
    
    ecc_vect = ecc_mat(1,:);
    
    lim_x = [min(ecc_vect)-0.05*max(ecc_vect-min(ecc_vect)),...
             max(ecc_vect)+0.05*max(ecc_vect-min(ecc_vect))];

    linewidth = 1;
    scFormat = 'r.';
    errFormat = '-k+';

    figure;
    latexify(40,18)

    subplot(1,3,1)
    hold on
    scatter(ecc_mat(:),rms(:),scFormat)
    errorbar(ecc_vect,mean(rms),std(rms),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Eccentricity, nd')
    ylabel('RMS Position Error, $km$')
    xlim(lim_x)
    setgrid

    subplot(1,3,2)
    hold on
    scatter(ecc_mat(:),period(:),scFormat)
    errorbar(ecc_vect,mean(period),std(period),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Eccentricity, nd')
    ylabel('Period Error, seconds')
    xlim(lim_x)
    setgrid

    subplot(1,3,3)
    hold on
    scatter(ecc_mat(:),ecc(:),scFormat)
    errorbar(ecc_vect,mean(ecc),std(ecc),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Eccentricity, nd')
    ylabel('Eccentricity Error, nd')
    xlim(lim_x)
    setgrid
    
    latexify(40,18,22)
    
    % clean up persistent variables
    clear rms period ecc ecc_mat

end

end
