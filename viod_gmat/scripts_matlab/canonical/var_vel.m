function var_vel(DU,TU,V,mu,noise,num_monte,pos_ref,period_ref,ecc_ref)
%VAR_VEL Plots error statistics for VIOD w.r.t mean velocity.
%
% Author:
%   Tiger Hou
%
 
% declare persistent variables
persistent rms period ecc v_mat

%% data collection mode
if nargin == 9
    
    noise = noise / DU * TU;
    
    rms_temp = nan(num_monte,1);
    period_temp = nan(num_monte,1);
    ecc_temp = nan(num_monte,1);
    
    % clean up unused velocity measurement slots
    V = vclean(V);
    vmean = mean(vecnorm(V,2,2));
    
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
    v_mat = [v_mat repmat(vmean,num_monte,1)];
    
    return;
    
%% data plotting mode
elseif nargin == 2
    
    v_vect = v_mat(1,:);
    
    lim_x = [min(v_vect)-0.05*max(v_vect-min(v_vect)),...
             max(v_vect)+0.05*max(v_vect-min(v_vect))];

    linewidth = 1;
    scFormat = 'r.';
    errFormat = '-k+';

    figure;
    latexify(40,18)

    subplot(1,3,1)
    hold on
    scatter(v_mat(:),rms(:),scFormat)
    errorbar(v_vect,mean(rms),std(rms),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Mean Velocity, DU/TU')
    ylabel('RMS Position Error, DU')
    xlim(lim_x)
    setgrid

    subplot(1,3,2)
    hold on
    scatter(v_mat(:),period(:),scFormat)
    errorbar(v_vect,mean(period),std(period),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Mean Velocity, DU/TU')
    ylabel('Period Error, TU')
    xlim(lim_x)
    setgrid

    subplot(1,3,3)
    hold on
    scatter(v_mat(:),ecc(:),scFormat)
    errorbar(v_vect,mean(ecc),std(ecc),errFormat,'LineWidth',linewidth)
    hold off
    xlabel('Mean Velocity, DU/TU')
    ylabel('Eccentricity Error, nd')
    xlim(lim_x)
    setgrid
    
    latexify(40,18,22)
    
    sgtitle(['DU = ' num2str(DU) ' $km$,  ',...
             'TU = ' num2str(TU) ' $s$,  ',...
             'DU/TU = ' num2str(DU/TU) ' km/s'])
    
    % clean up persistent variables
    clear rms period ecc v_mat

end

end
