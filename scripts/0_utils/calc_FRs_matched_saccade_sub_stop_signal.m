function [event_frs,sdfs,eventRasters] = calc_FRs_matched_saccade_sub_stop_signal(spikeInds,...
    baseline_win,event_alignment,matchedTrls,window,time,step,subSDF,RT,SS,SSRT)
% Calculates firing rates and spike density functions for matched trials
% but subtracts the stop signal response from the SSRT aligned cancelled vs
% no stop
% Inputs: spikeInds - spike timing data
%         baseline_win - window for baseline firing rate
%         event_alignment - event to align spikes to
%         ssd_offsets - stop signal delays
%         matchedTrls - matched trial indices
%         window/time/step - parameters for calculation windows
%         subSDF - the difference between cancelled and no stop trials
%           aligned to SSRT
%         RT - per trial reaction times
%         SS - per trial stop-signal times
%         SSRT - the average SSRT

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
sdfs = nan(size(rastersSac,1),2001);
eventRasters = nan(size(rastersSac,1),2001);
tempRasters = nan(size(rastersSac,1),3001);
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
           eventRasters(currTrl,:) = rastersSac(currTrl,1:1 + 2000);
           sdfs(currTrl,:) = mL_SDF(eventRasters(currTrl,1:2001));
           
           % Calculate Z-scored firing rates
           for fr_i = 1:length(steps)
               event_frs(currTrl,fr_i) = (sum(eventRasters(1,steps(fr_i):steps(fr_i)+time))...
                   ./(time/1000)-mean_fr)./std_fr;
           end
       end
       
       % Process non-cancelled trials
       for mtrl = 1:length(matchedTrls{ssd_i,2})
           currTrl = matchedTrls{ssd_i,2}(mtrl);
           %calculate the SSRT aligned trial SDF
           tempRasters(currTrl,:) = rastersTar(currTrl,1:1 + 3000);
           ssrtSDF = mL_SDF(tempRasters(currTrl,1:3001));
           %subtract the canc-nostop difference (still aligned to target
           correctedSDF = ssrtSDF - [zeros(1,SSRT + SS(currTrl)) subSDF zeros(1,1000-(SSRT + SS(currTrl)))];
           %realign to actual RT
           sdfs(currTrl,:) = correctedSDF(RT(currTrl):RT(currTrl)+2000);
           
           % Calculate Z-scored firing rates
           for fr_i = 1:length(steps)
               event_frs(currTrl,fr_i) = (sum(tempRasters(1,steps(fr_i):steps(fr_i)+time))...
                   ./(time/1000)-mean_fr)./std_fr;
           end
       end
   end
end