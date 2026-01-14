function [event_frs,sdfs] = calc_FRs_and_zscore(spikeInds,baseline_win,event_alignment,...
    window,time,step)

steps =( window(1):step:window(2)-time)+1000;
rastersTar = generateRasterDaJo(spikeInds.target);
mean_fr = nanmean(sum(rastersTar(:,1000+baseline_win),2))/(length(baseline_win)/1000);
std_fr = std(sum(rastersTar(:,1000+baseline_win),2)/(length(baseline_win)/1000));

rastersEvent = generateRasterDaJo(spikeInds.(event_alignment));
event_frs = zeros(size(rastersEvent,1),length(steps));
for trial_i = 1:size(spikeInds.target,1)
    for fr_i = 1:length(steps)
    event_frs(trial_i,fr_i) = (sum(rastersEvent(trial_i,steps(fr_i):...
        steps(fr_i)+time))./(time/1000)-mean_fr)./std_fr;
    end
end
sdfs = cell2mat(arrayfun(@(x) mL_SDF(rastersEvent(x,:)),1:size(rastersEvent,1),'uni',0)');
end