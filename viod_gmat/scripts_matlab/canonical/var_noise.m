function var_noise(V,mu,DU,TU,num,pos_ref,period_ref,ecc_ref)
%VAR_NOISE Plots error statistics for VIOD w.r.t noise stdev.
%
% Author:
%   Tiger Hou
%

% define noise
noise = [0.01 0.1 0.5 1 2 3 5 7 10] / DU * TU;

% clean up unused velocity measurement slots
V = vclean(V);

% RMS position error
rms = nan(num,length(noise));
% orbit period error
period = nan(num,length(noise));
% eccentricity error
ecc = nan(num,length(noise));

for i = 1:length(noise)
    
    for j = 1:num
        
        vn = addnoise(V,noise(i));
        rn = viod(vn,mu);
        rr = rn(1,:);
        vv = vn(1,:);
        
        rms(j,i)    = norm(rr-pos_ref);
        sma_est     = norm(rr)/(2 - norm(rr)*dot(vv,vv)/mu);
        period(j,i) = 2*pi*sqrt(sma_est^3/mu) - period_ref;
        ecc(j,i)    = norm(((dot(vv,vv)-mu/norm(rr))*rr-dot(rr,vv)*vv)/mu)...
                    - ecc_ref;
        
    end
    
end

noise = noise / 1000;
noise_mat = repmat(noise,num,1);

%% plotting

linewidth = 1;
scFormat = 'r.';
errFormat = '-k+';

lim_x = [min(noise)-0.05*max(noise-min(noise)),...
         max(noise)+0.05*max(noise-min(noise))];

figure;
latexify(40,18,22)

subplot(1,3,1)
hold on
scatter(noise_mat(:),rms(:),scFormat)
errorbar(noise,mean(rms),std(rms),errFormat,'LineWidth',linewidth)
hold off
xlabel('Noise, DU/TU')
ylabel('Position Error, DU')
xlim(lim_x)
setgrid

subplot(1,3,2)
hold on
scatter(noise_mat(:),period(:),scFormat)
errorbar(noise,mean(period),std(period),errFormat,'LineWidth',linewidth)
hold off
xlabel('Noise, DU/TU')
ylabel('Period Error, TU')
xlim(lim_x)
setgrid

subplot(1,3,3)
hold on
scatter(noise_mat(:),ecc(:),scFormat)
errorbar(noise,mean(ecc),std(ecc),errFormat,'LineWidth',linewidth)
hold off
xlabel('Noise, DU/TU')
ylabel('Eccentricity Error, nd')
xlim(lim_x)
setgrid
sgtitle(['DU = ' num2str(DU) ' $km$,  ',...
         'TU = ' num2str(TU) ' $s$,  ',...
         'DU/TU = ' num2str(DU/TU) ' km/s'])
    
latexify(40,18,22)

end
