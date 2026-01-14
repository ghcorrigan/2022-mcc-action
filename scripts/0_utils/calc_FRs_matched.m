function [event_frs,sdfs,eventRasters] = calc_FRs_matched(spikeInds,...
    baseline_win,event_alignment,ssd_offsets,matchedTrls,window,time,step)
% Calculates firing rates and spike density functions for matched trials
% Inputs: spikeInds - spike timing data
%         baseline_win - window for baseline firing rate
%         event_alignment - event to align spikes to
%         ssd_offsets - stop signal delays
%         matchedTrls - matched trial indices
%         window/time/step - parameters for calculation windows

% Calculate time steps for analysis windows
steps = (window(1):step:window(2)-time)+1000;

% Generate target-aligned spike rasters
rastersTar = generateRasterDaJo(spikeInds.target);

% Calculate baseline firing rate stats
mean_fr = nanmean(sum(rastersTar(:,1000+baseline_win),2))/(length(baseline_win)/1000);
std_fr = std(sum(rastersTar(:,1000+baseline_win),2)/(length(baseline_win)/1000));

% Initialize output arrays
sdfs = nan(size(rastersTar,1),2001);
eventRasters = nan(size(rastersTar,1),2001);
rastersEvent = generateRasterDaJo(spikeInds.(event_alignment));
event_frs = zeros(size(rastersEvent,1),length(steps));

% Calculate event-aligned firing rates
for trial_i = 1:size(spikeInds.target,1)
   for fr_i = 1:length(steps)
       event_frs(trial_i,fr_i) = (sum(rastersEvent(trial_i,steps(fr_i):steps(fr_i)+time))...
           ./(time/1000)-mean_fr)./std_fr;
   end
end

% Process matched trials
for ssd_i = 1:length(matchedTrls)
   if ~isempty(matchedTrls{ssd_i,1})
       % Process non-cancelled trials
       for mtrl = 1:length(matchedTrls{ssd_i,1})
           currTrl = matchedTrls{ssd_i,1}(mtrl);
           eventRasters(currTrl,:) = rastersTar(currTrl,ssd_offsets(ssd_i)+1:ssd_offsets(ssd_i)+1 + 2000);
           sdfs(currTrl,:) = mL_SDF(eventRasters(currTrl,1:2001));
           
           % Calculate Z-scored firing rates
           for fr_i = 1:length(steps)
               event_frs(currTrl,fr_i) = (sum(eventRasters(1,steps(fr_i):steps(fr_i)+time))...
                   ./(time/1000)-mean_fr)./std_fr;
           end
       end
       
       % Process cancelled trials 
       for mtrl = 1:length(matchedTrls{ssd_i,2})
           currTrl = matchedTrls{ssd_i,2}(mtrl);
          eventRasters(currTrl,:) = rastersTar(currTrl,ssd_offsets(ssd_i)+1:ssd_offsets(ssd_i)+1 + 2000);
           sdfs(currTrl,:) = mL_SDF(eventRasters(currTrl,1:2001));
           
           % Calculate Z-scored firing rates
           for fr_i = 1:length(steps)
               event_frs(currTrl,fr_i) = (sum(eventRasters(1,steps(fr_i):steps(fr_i)+time))...
                   ./(time/1000)-mean_fr)./std_fr;
           end
       end
   end
end