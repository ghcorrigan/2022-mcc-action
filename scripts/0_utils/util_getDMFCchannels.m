function dfmc_ch_mapping = util_getDMFCchannels(logInfo,session_i)

ctx_ch = logInfo.ctx_top_ch(session_i);
if isnan(ctx_ch)
   ctx_ch = NaN; 
end


acc_ch(:,1) = 1:32;
acc_ch(:,2) = zeros(length(1:32),1);

walk = 0;
while ctx_ch-(walk+1) > 0
    walk = walk + 1;
    acc_ch(ctx_ch-(walk),2) = -walk;
end
walk = 0;
while ctx_ch+(walk+1) <= 32 
    walk = walk + 1;
    acc_ch(ctx_ch+(walk),2) = walk;
end

acc_ch(:,3) = acc_ch(:,2) * logInfo.electrode_spc(session_i);

n_dMCC_ch = sum(acc_ch(:,2) < 0);
n_vMCC_ch = sum(acc_ch(:,2) > 0);

if n_dMCC_ch > 0 & n_vMCC_ch > 0
labels = [repmat({'notCtx'},n_dMCC_ch,1);...
    repmat({'DMFC'},n_vMCC_ch+1,1)];
elseif n_dMCC_ch == 0 & n_vMCC_ch > 0
    labels = [{'Ctx'};repmat({'DMFC'},n_vMCC_ch,1)];

else
    labels = repmat({'?'},32,1);
end


dfmc_ch_mapping = table(acc_ch(:,1), acc_ch(:,2), acc_ch(:,3), labels,...
    'VariableNames',{'channel','ch_depth_acs','um_depth_acs','area'});



end
