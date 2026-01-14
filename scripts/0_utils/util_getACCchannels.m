function acs_ch_mapping = util_getACCchannels(logInfo,session_i)

acs_ch = logInfo.acs_ch(session_i);
% if acs_ch == 0
%    acs_ch = NaN; 
% end


acc_ch(:,1) = 1:32;
acc_ch(:,2) = zeros(length(1:32),1);

walk = 0;
while acs_ch-(walk+1) > 0
    walk = walk + 1;
    acc_ch(acs_ch-(walk),2) = -walk;
end
walk = 0;
while acs_ch+(walk+1) <= 32 
    walk = walk + 1;
    acc_ch(acs_ch+(walk),2) = walk;
end

acc_ch(:,3) = acc_ch(:,2) * logInfo.electrode_spc(session_i);

n_dMCC_ch = sum(acc_ch(:,2) < 0);
n_vMCC_ch = sum(acc_ch(:,2) > 0);

if n_dMCC_ch > 0 & n_vMCC_ch > 0
labels = [repmat({'dMCC'},n_dMCC_ch,1);...
    {'ACS'};repmat({'vMCC'},n_vMCC_ch,1)];
elseif n_dMCC_ch == 0 & n_vMCC_ch > 0
    labels = [{'ACS'};repmat({'vMCC'},n_vMCC_ch,1)];
    if length(labels)>32
        labels = labels(end-31:end);
    end
elseif n_vMCC_ch == 0 & n_dMCC_ch > 0
    labels = [repmat({'dMCC'},n_dMCC_ch,1);{'ACS'}];
    if length(labels)>32
        labels(33:end)=[];
    end
else
    labels = repmat({'?'},32,1);
end


acs_ch_mapping = table(acc_ch(:,1), acc_ch(:,2), acc_ch(:,3), labels,...
    'VariableNames',{'channel','ch_depth_acs','um_depth_acs','area'});



end
