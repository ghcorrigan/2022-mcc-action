function [event_frs,mean_preTar] = calc_FRs_no_zscore(spikeInds,baseline_win,event_alignment,...
    window,time,step)

steps =( window(1):step:window(2)-time)+1000;
rastersTar = generateRasterDaJo(spikeInds.target);
mean_preTar = nanmean(sum(rastersTar(:,1000+baseline_win),2))/(length(baseline_win)/1000);

rastersEvent = generateRasterDaJo(spikeInds.(event_alignment));
event_frs = zeros(size(rastersEvent,1),length(steps));
for trial_i = 1:size(spikeInds.target,1)
    for fr_i = 1:length(steps)
    event_frs(trial_i,fr_i) = (sum(rastersEvent(trial_i,steps(fr_i):...
        steps(fr_i)+time))./(time/1000));
    end
end
end