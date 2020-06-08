function var_obsv(V,mu,noise,num_monte,num_obsv,pos_ref,period_ref,ecc_ref)
%VAR_OBSV Plots error statistics for VIOD w.r.t number of observations.
%
% Author:
%   Tiger Hou
%

% declare persistent variables
persistent rms period ecc obsv_mat

%% data collection mode
if nargin == 8
    
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
    obsv_mat = [obsv_mat repmat(num_obsv,num_monte,1)];
    
    return;
    
%% data plotting mode
else
    
    obsv = obsv_mat(1,:);
    
    lim_x = [min(obsv)-0.05*max(obsv-min(obsv)),...
             max(obsv)+0.05*max(obsv-min(obsv))];

    linewidth = 1;
    scFormat = 'r.';
    errFormat = '-k+';

    figure;
    latexify(40,18)

    subplot(1,3,1)
    hold on
    scatter(obsv_mat(:),rms(:),scFormat)
    errorbar(obsv,mean(rms),std(rms),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Number of Observations, nd')
    ylabel('RMS Position Error, $km$')
    xlim(lim_x)
    setgrid

    subplot(1,3,2)
    hold on
    scatter(obsv_mat(:),period(:),scFormat)
    errorbar(obsv,mean(period),std(period),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Number of Observations, nd')
    ylabel('Period Error, seconds')
    xlim(lim_x)
    setgrid

    subplot(1,3,3)
    hold on
    scatter(obsv_mat(:),ecc(:),scFormat)
    errorbar(obsv,mean(ecc),std(ecc),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Number of Observations, nd')
    ylabel('Eccentricity Error, nd')
    xlim(lim_x)
    setgrid
    
    latexify(40,18,22)
    
    % clean up persistent variables
    clear rms period ecc obsv_mat

end

end
