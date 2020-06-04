function var_vel(V,mu,noise,num_monte,pos_ref,period_ref,ecc_ref)
%VAR_VEL Plots error statistics for VIOD w.r.t mean velocity.
%
% Author:
%   Tiger Hou
%
 
% declare persistent variables
persistent rms period ecc v_mat

%% data collection mode
if nargin == 7
    
    rms_temp = nan(num_monte,1);
    period_temp = nan(num_monte,1);
    ecc_temp = nan(num_monte,1);
    
    % clean up unused velocity measurement slots
    V = vclean(V);
    vmean = mean(vecnorm(V,2,2));
    
    for i = 1:num_monte

        vn = addnoise(V,noise);
        rn = viod(vn,mu);
        rr = rn(end,:);
        vv = vn(end,:);

        rms_temp(i)    = norm(pos_ref-rr);
        sma_est        = norm(rr)/(2 - norm(rr)*dot(vv,vv)/mu);
        period_temp(i) = 2*pi*sqrt(sma_est^3/mu) - period_ref;
        ecc_temp(i)    = norm(((dot(vv,vv)-mu/norm(rr))*rr-dot(rr,vv)*vv)/mu)...
                       - ecc_ref;

    end
    
    rms = [rms rms_temp];
    period = [period period_temp];
    ecc = [ecc ecc_temp];
    v_mat = [v_mat repmat(vmean,num_monte,1)];
    
    return;
    
%% data plotting mode
else
    
    v_vect = v_mat(1,:);

    linewidth = 1;
    scFormat = 'r.';
    errFormat = '-k+';

    figure;
    latexify(45,18)

    subplot(1,3,1)
    hold on
    scatter(v_mat(:),rms(:),scFormat)
    errorbar(v_vect,mean(rms),std(rms),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Mean Velocity, $km/s$')
    ylabel('RMS Position Error, $km$')

    subplot(1,3,2)
    hold on
    scatter(v_mat(:),period(:),scFormat)
    errorbar(v_vect,mean(period),std(period),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Mean Velocity, $km/s$')
    ylabel('Period Error, seconds')

    subplot(1,3,3)
    hold on
    scatter(v_mat(:),ecc(:),scFormat)
    errorbar(v_vect,mean(ecc),std(ecc),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Mean Velocity, $km/s$')
    ylabel('Eccentricity Error, nd')
    
    % clean up persistent variables
    clear rms period ecc v_mat

end

end
