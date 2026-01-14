function [event_frs,sdfs,eventRasters] = calc_FRs_matched_saccade(spikeInds,...
    baseline_win,event_alignment,matchedTrls,window,time,step)
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

% Generate saccade-aligned spike rasters
rastersSac = generateRasterDaJo(spikeInds.saccade);

% Calculate baseline firing rate stats
mean_fr = nanmean(sum(rastersTar(:,1000+baseline_win),2))/(length(baseline_win)/1000);
std_fr = std(sum(rastersTar(:,1000+baseline_win),2)/(length(baseline_win)/1000));

% Initialize output arrays
sdfs = nan(size(rastersSac,1),3001);
eventRasters = nan(size(rastersSac,1),3001);
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
for ssd_i = 1:size(matchedTrls,1)
   if ~isempty(matchedTrls{ssd_i,1})
       % Process no stop
       for mtrl = 1:length(matchedTrls{ssd_i,1})
           currTrl = matchedTrls{ssd_i,1}(mtrl);
           eventRasters(currTrl,:) = rastersSac(currTrl,1:1 + 3000);
           sdfs(currTrl,:) = mL_SDF(eventRasters(currTrl,1:3001));
           
           % Calculate Z-scored firing rates
           for fr_i = 1:length(steps)
               event_frs(currTrl,fr_i) = (sum(eventRasters(1,steps(fr_i):steps(fr_i)+time))...
                   ./(time/1000)-mean_fr)./std_fr;
           end
       end
       
       % Process non-cancelled trials
       for mtrl = 1:length(matchedTrls{ssd_i,2})
           currTrl = matchedTrls{ssd_i,2}(mtrl);
          eventRasters(currTrl,:) = rastersSac(currTrl,1:1 + 3000);
           sdfs(currTrl,:) = mL_SDF(eventRasters(currTrl,1:3001));
           
           % Calculate Z-scored firing rates
           for fr_i = 1:length(steps)
               event_frs(currTrl,fr_i) = (sum(eventRasters(1,steps(fr_i):steps(fr_i)+time))...
                   ./(time/1000)-mean_fr)./std_fr;
           end
       end
   end
end