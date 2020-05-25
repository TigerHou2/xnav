function [data_ell,data_hyp,mean_ell,std_ell,...
          mean_hyp,std_hyp] = Data_Consolidate_Additive(data,init_cond)
%DATA_CONSOLIDATE Summary of this function goes here
%   Detailed explanation goes here

obsv_min = init_cond.min;
obsv_max = init_cond.max;
trials   = init_cond.num;

data_ell = [];
data_hyp = [];
mean_ell = [];
mean_hyp = [];
std_ell  = [];
std_hyp  = [];
temp_ell = [];
temp_hyp = [];
prev_int = obsv_min;

for i = 1:((length(obsv_min:obsv_max))*trials)
    if not(data(i).int == prev_int)
        if size(temp_ell,1) <= 1
            mean_ell = [mean_ell; temp_ell];
            std_ell  = [std_ell ; temp_ell];
        else
            mean_ell = [mean_ell; mean(temp_ell)];
            std_ell  = [std_ell ; std(temp_ell)];
        end
        temp_ell = [];
        if size(temp_hyp,1) <= 1
            mean_hyp = [mean_hyp; temp_hyp];
            std_hyp  = [std_hyp ; temp_hyp];
        else
            mean_hyp = [mean_hyp; mean(temp_hyp)];
            std_hyp  = [std_hyp ; std(temp_hyp)];
        end
        temp_hyp = [];
        prev_int = data(i).int;
    end
    if data(i).is_elliptic
        temp_ell = [temp_ell; [data(i).int data(i).df data(i).orb]];
        data_ell = [data_ell; [data(i).int data(i).df data(i).orb]];
    else
        temp_hyp = [temp_hyp; [data(i).int data(i).df data(i).orb]];
        data_hyp = [data_hyp; [data(i).int data(i).df data(i).orb]];
    end
end

end

