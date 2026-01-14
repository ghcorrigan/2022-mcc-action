function [glm_output, encoding_flag, encoding_beta,inh_pnc,frs,sig_bins,pvals] = ...
    stopping_glm_func_fr_aligned(neuron_i, mcc_map_info, behData)
rng(1988)
neuralFilename = mcc_map_info.session{neuron_i};
neuronLabel = mcc_map_info.unit{neuron_i};
% fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))


% Load in pre-processed spike data
% import_data = struct(); import_data = load_behFile(dirs,behFilename);
regAlpha = .02;
doPerm = 0;
load(['E:\conflict\results\DATA_SPK_STATIONARITY\Ben\stationarityData_' mcc_map_info.probeName{neuron_i}])
load(['E:\conflict\results\probeTrls\neuronData_' mcc_map_info.probeName{neuron_i}])
tempInd = find(contains(neuronData.DSPname,neuronLabel));
stationInd = find(stationaritySelection.included_trials.neuronIdx == tempInd);

%% Extract: get relevant data for GLM table

reg_tbl = table;

reg_tbl.trial_type = zeros(length(behData.trialEventTimes{1}.isNoStop),1);
reg_tbl.trial_type(behData.trialEventTimes{1}.isNoStop==1) = 1;
reg_tbl.trial_type(behData.trialEventTimes{1}.isCanceled==1) = 2;
reg_tbl.trial_type(behData.trialEventTimes{1}.isNonCanceled==1) = 3;
reg_tbl.seen(~isnan(behData.trialEventTimes{1}.stopSignal)) = 1;
reg_tbl.ssd = behData.trialEventTimes{1}.ssdInd+1;
reg_tbl.trial_number = [1:height(reg_tbl)]';
% for trial_i = 1:size(import_data.events.stateFlags_,1)
%     
%     % Trial type ---------------------------
%     % - No-stop trials
%     if import_data.events.stateFlags_.IsGoCorrect(trial_i) == 1
%         reg_tbl.trial_type(trial_i) = 1;
% %         reg_tbl.RT(trail_i) = import_data.events.stateFlags_.
%     % - Canceled trials
%     elseif import_data.events.stateFlags_.IsCancel(trial_i) == 1
%         reg_tbl.trial_type(trial_i) = 2;
%     % - Non-canceled trials
%     elseif import_data.events.stateFlags_.IsNonCancelledNoBrk(trial_i) == 1 ||...
%             import_data.events.stateFlags_.IsNonCancelledBrk(trial_i) == 1
%         reg_tbl.trial_type(trial_i) = 3;
%         if isnan(behData.trialEventTimes{1}.stopSignal(trial_i))
%             %unseen
%             reg_tbl.seen(trial_i) = 0;
%         else
%             %seen
%             reg_tbl.seen(trial_i) = 1;
%             
%         end
%     % - Abort trials
%     else
%         reg_tbl.trial_type(trial_i) = 0;
%     end
%     
%     % Stop signal delay ---------------------------
%     reg_tbl.ssd(trial_i) = import_data.events.stateFlags_.UseSsdIdx(trial_i)+1;
%     
%     % Value ---------------------------------------
%     if import_data.events.stateFlags_.IsLoRwrd(trial_i) == 1
%         reg_tbl.value(trial_i) = 1;
%     else
%         reg_tbl.value(trial_i) = 2;
%     end
%     
%     % Trial Number ---------------------------------------
%     reg_tbl.trial_number(trial_i) = trial_i;
%     
% end
goodtrls = false(height(reg_tbl),1);
goodtrls(stationaritySelection.included_trials.trials{stationInd}) = true;
%check what the shortest SSD is, and if shorter than 80ms, it gets removed
if behData.stopSignalBeh.inh_SSD(1)<80
goodtrls(reg_tbl.ssd==1 ) = false;
end
% tooFew = find(cellfun(@length, behData.stopSignalBeh.ssd_ttx.C)<12);
% if ~isempty(tooFew)
% for ii = 1:length(tooFew)
%     goodtrls(reg_tbl.ssd==tooFew(ii)-1) = false;
% end
% end
reg_tbl.pnc = zeros(height(reg_tbl),1);
for ii = 1:max(reg_tbl.ssd)
    reg_tbl.pnc(reg_tbl.ssd==ii) = behData.stopSignalBeh.inh_pnc(ii);
    inh_pnc(1,ii) = behData.stopSignalBeh.inh_pnc(ii);
    inh_pnc(2,ii) = behData.stopSignalBeh.inh_SSD(ii)';
end
%intialize fr struct
frs.canc = {};
    frs.ssd = {};
    frs.pnc = {};
if height(reg_tbl(goodtrls,:))==0 || sum(reg_tbl.trial_type(goodtrls)==2)<30
    glm_output.trial_type.sig_times = [];
    glm_output.trial_type.beta_weights = [];
    glm_output.trial_type.permPos =[];
    glm_output.trial_type.pvals =[];
    glm_output.ssd.sig_times = [];
    glm_output.ssd.sig_timePos = [];
    glm_output.ssd.beta_weights = [];
    glm_output.ssd.permPos = [];
    glm_output.ssd.permSig = [];
    glm_output.ssd.pvals = [];
    glm_output.seen.sig_times = [];
    glm_output.seen.beta_weights = [];
    glm_output.pnc.sig_times = [];
    glm_output.pnc.sig_timePos = [];
    glm_output.pnc.beta_weights = [];
    glm_output.pnc.permPos =[];
    glm_output.pnc.pvals =[];
    glm_output.nonCanc.sig_times = [];
    glm_output.nonCanc.beta_weights = [];
    encoding_flag = zeros(1,6); encoding_beta = zeros(1,5);
    inh_pnc = {0};
    sig_bins = cell(1,3);
    pvals = struct();
    remove = 1;
    return
end
%% Setup spike data into GLM
spk_data_in = load(fullfile('E:\conflict\data\','Spikes',...
    [neuralFilename '_Spikes_' neuronLabel '.mat']));
window_shift = 10;
fr_window = 100;
event_alignment = 'ssrt';
baseline_win = [-600:-500];
[event_frs,sdfs] = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
    [-50, 350],fr_window,window_shift);

[~, n_times] = size(event_frs); % Get the number of windows in the time averaged data

reg_tbl.win_fr = event_frs; % Add the firing rate over this whole window to the GLM table

%% Run GLM: trial type
glm_output.trial_type.sig_times = []; 
glm_output.trial_type.beta_weights = []; 
glm_output.trial_type.permPos = zeros(1,n_times); 
glm_output.trial_type.pvals = zeros(1,n_times); 
u_t_mdl = [];
reg_tbl.is_canc = reg_tbl.trial_type==2;
reg_tbl.trl_numZ = zscore(reg_tbl.trial_number);
%match trials to calculate noStop firing rates
ssrt_delays = round(behData.stopSignalBeh.inh_SSD + ...
    behData.stopSignalBeh.ssrt.integrationWeighted);

trlInds = reg_tbl.trial_number;
goodTrlNums = find(goodtrls);
[matchedInds,matchedBySSD] = trlMatcherCanc(behData,goodTrlNums);
event_frs_matched = calc_FRs_matched(spk_data_in.Spikes,baseline_win,event_alignment,...
    ssrt_delays,matchedBySSD,[-50, 350],fr_window,window_shift);
reg_tbl.win_fr = event_frs_matched;
%assign the matched SSD
for ssd_i = 1:length(matchedBySSD)
    reg_tbl.ssd(matchedBySSD{ssd_i,1}) = ssd_i;
end
matchedInds = trlInds(ismember(trlInds,matchedInds));

reg_tbl_trialtype = reg_tbl(matchedInds,:);
reg_tbl_nonCanc = reg_tbl(reg_tbl.trial_type == 2 | reg_tbl.trial_type == 3,:);
ttypePerm = zeros(height(reg_tbl_trialtype),500);
for ii = 1:500
    ttypePerm(:,ii) = randsample(reg_tbl_trialtype.is_canc, height(reg_tbl_trialtype));
end
% For each averaged time point
for timepoint_i = 1:n_times
    
    % Input the timepoint specific firing times
    reg_tbl_trialtype.firing_rate =  reg_tbl_trialtype.win_fr(:,timepoint_i);
    
    % Fit the GLM for this
   
    u_t_mdl = fitlm(reg_tbl_trialtype,'firing_rate ~ is_canc + trl_numZ');
    if doPerm ==1
    for jj = 1:500
        reg_tbl_trialtype.isCancP = ttypePerm(:,jj);
        perm_t_mdl = fitlm(reg_tbl_trialtype,'firing_rate ~ isCancP + trl_numZ');
        coeffs(jj) = perm_t_mdl.Coefficients.Estimate(3);
    end
   
    permLoc = length(find(coeffs<u_t_mdl.Coefficients.Estimate(2)));
        glm_output.trial_type.permPos(1,timepoint_i) = permLoc; % trial type

    end
    % GLM output ---------------------------------
    % - Trial type
    glm_output.trial_type.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % trial type
    glm_output.trial_type.pvals(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2); % trial type
    glm_output.trial_type.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.Estimate(2); % trial type

    % - Trial number
    glm_output.trial_type.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.trial_type.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.Estimate(3); % trial_number
    
end

%% Run GLM: stop_signal delay
%have to subtract the average activity for atency matched no stop
for ssd_i = 1:length(matchedBySSD)
   noStopAvg(ssd_i,:) = mean( reg_tbl_trialtype.win_fr(...
       reg_tbl_trialtype.ssd == ssd_i & reg_tbl_trialtype.trial_type ==1,:));
   reg_tbl_trialtype.fr_diff(reg_tbl_trialtype.ssd==ssd_i,:) = ...
       reg_tbl_trialtype.win_fr(reg_tbl_trialtype.ssd==ssd_i,:) -...
       noStopAvg(ssd_i,:);
   reg_tbl_trialtype.ssdVal(reg_tbl_trialtype.ssd==ssd_i,:) = inh_pnc(2,ssd_i);
   reg_tbl_trialtype.pnc(reg_tbl_trialtype.ssd==ssd_i,:) = inh_pnc(2,ssd_i); 
end
glm_output.ssd.sig_times = zeros(1,n_times); 
glm_output.ssd.beta_weights = zeros(1,n_times); 
glm_output.ssd.permPos = zeros(1,n_times);
glm_output.ssd.permSig = zeros(1,n_times);
glm_output.ssd.pvals = zeros(1,n_times);
u_t_mdl = [];

reg_tbl_stopping = reg_tbl_trialtype(reg_tbl_trialtype.is_canc==1,:);
%check that there are enough trials of eac SSD
usedSSDs = unique(reg_tbl_stopping.ssd);
for ssd_i = 1:length(usedSSDs)
    if height(reg_tbl_stopping(reg_tbl_stopping.ssd==usedSSDs(ssd_i),:))<5
        reg_tbl_stopping(reg_tbl_stopping.ssd==usedSSDs(ssd_i),:) = [];
    end
end
if doPerm==1
sigPermTimes = find(abs(glm_output.trial_type.permPos-250)>238 |...
    glm_output.trial_type.pvals < regAlpha);
else
    sigPermTimes = find(glm_output.trial_type.pvals < regAlpha);
end
ssdPerm = zeros(height(reg_tbl_stopping),500);
for ii = 1:500
    ssdPerm(:,ii) = randsample(reg_tbl_stopping.ssdVal, height(reg_tbl_stopping));
end
% For each averaged time point
for timepoint_i = sigPermTimes
    
    % Input the timepoint specific firing times
    reg_tbl_stopping.firing_rate = reg_tbl_stopping.fr_diff(:,timepoint_i); 
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ ssdVal + trl_numZ');
    if doPerm ==1
    for jj = 1:500
        reg_tbl_stopping.ssdP = ssdPerm(:,jj);
        perm_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ ssdP + trl_numZ');
        coeffs(jj) = perm_t_mdl.Coefficients.Estimate(3);
        glm_output.ssd.permPos(1,timepoint_i) = permLoc;
        glm_output.ssd.permSig(1,timepoint_i) = abs(permLoc-250)>=238;

    end
   
    permLoc = length(find(coeffs<u_t_mdl.Coefficients.Estimate(2)));
    end
    % GLM output ---------------------------------   
    % - Stop-signal delay
    glm_output.ssd.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < regAlpha; % ssd
    glm_output.ssd.pvals(1,timepoint_i) = u_t_mdl.Coefficients.pValue(3); % ssd
    glm_output.ssd.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.Estimate(3); % ssd
    
%     % - Value
%     glm_output.ssd.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % value
%     glm_output.ssd.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % value
%         
%     % - SSD x Value
%     glm_output.ssd.sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(5) < .01; % ssd x value
%     glm_output.ssd.beta_weights(3,timepoint_i) = u_t_mdl.Coefficients.tStat(5); % % ssd x value
    
    % - Trial number
    glm_output.ssd.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.ssd.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.Estimate(3); % trial_number

    
end
%% Run GLM: pNC
glm_output.pnc.sig_times = zeros(1,n_times); 
glm_output.pnc.beta_weights = zeros(1,n_times); 
glm_output.pnc.permPos = zeros(1,n_times);
glm_output.pnc.permSig = zeros(1,n_times);
glm_output.pnc.pvals = zeros(1,n_times);
u_t_mdl = [];

% if not doing perm testing, only test pvals
if doPerm
sigPermTimes = find(abs(glm_output.trial_type.permPos-250)>238 |...
    glm_output.trial_type.pvals < regAlpha);
else
    sigPermTimes = find(glm_output.trial_type.pvals < regAlpha);
end
ssdPerm = zeros(height(reg_tbl_stopping),500);
for ii = 1:500
    ssdPerm(:,ii) = randsample(reg_tbl_stopping.ssd, height(reg_tbl_stopping));
end
% For each averaged time point
for timepoint_i = sigPermTimes
    
    % Input the timepoint specific firing times
    reg_tbl_stopping.firing_rate = reg_tbl_stopping.fr_diff(:,timepoint_i); 
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ pnc + trl_numZ');
    if doPerm ==1
    for jj = 1:500
        reg_tbl_stopping.ssdP = ssdPerm(:,jj);
        perm_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ ssdP + trl_numZ');
        coeffs(jj) = perm_t_mdl.Coefficients.Estimate(3);
    end
   
    permLoc = length(find(coeffs<u_t_mdl.Coefficients.Estimate(2)));
    glm_output.pnc.permPos(1,timepoint_i) = permLoc;
    glm_output.pnc.permSig(1,timepoint_i) = abs(permLoc-250)>=238;
    end
    % GLM output ---------------------------------   
    % - Stop-signal delay
    glm_output.pnc.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % ssd
    glm_output.pnc.pvals(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2); % ssd
    glm_output.pnc.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.Estimate(2); % ssd
    
%     % - Value
%     glm_output.ssd.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % value
%     glm_output.ssd.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % value
%         
%     % - SSD x Value
%     glm_output.ssd.sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(5) < .01; % ssd x value
%     glm_output.ssd.beta_weights(3,timepoint_i) = u_t_mdl.Coefficients.tStat(5); % % ssd x value
    
    % - Trial number
    glm_output.pnc.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.pnc.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.Estimate(3); % trial_number

    
end
%% determine if there is a difference between cancelled and noncancelled
glm_output.nonCanc.sig_times = []; 
glm_output.nonCanc.beta_weights = []; 
u_t_mdl = [];
matchedInds = trlMatcherCancNC(behData,goodTrlNums);

reg_tbl_canc = reg_tbl(matchedInds,:);

for timepoint_i = 1:n_times
    
    % Input the timepoint specific firing times
    reg_tbl_canc.firing_rate = reg_tbl_canc.win_fr(:,timepoint_i);    
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_canc,'firing_rate ~ is_canc + trl_numZ');
    
    % GLM output ---------------------------------
    % - Trial type
    glm_output.nonCanc.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % trial type
    glm_output.nonCanc.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % trial type
    
    % - Trial number
    glm_output.nonCanc.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.nonCanc.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % trial_number
    
end
%% seen vs unseen
glm_output.seen.sig_times = []; 
glm_output.seen.beta_weights = []; 
u_t_mdl = []; 
trlInds = reg_tbl.trial_number;
matchedInds = trlMatcherNC(behData,goodTrlNums);
matchedInds = trlInds(ismember(trlInds,matchedInds));

reg_tbl_seen = reg_tbl(matchedInds,:);

for timepoint_i = 1:n_times
    
    % Input the timepoint specific firing times
    reg_tbl_seen.firing_rate = reg_tbl_seen.win_fr(:,timepoint_i);  
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_seen,'firing_rate ~ seen + trl_numZ');
    
    % GLM output ---------------------------------
    % - Trial type
    glm_output.seen.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % trial type
    glm_output.seen.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % trial type
    
    % - Trial number
    glm_output.seen.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.seen.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % trial_number
    
end

%% Determine periods of significance
%if seen isn't significant, we don't care about SSD being significant
% glm_output.ssd.sig_times(1,glm_output.seen.sig_times(1,:) ==0)=0;
% 
%if conflict isn't significant, we don't care about SSD being significant
glm_output.ssd.sig_times(1,glm_output.trial_type.sig_times(1,:) ==0)=0;
glm_output.ssd.sig_timePos = glm_output.ssd.sig_times & (glm_output.ssd.beta_weights>0);
glm_output.pnc.sig_times(1,glm_output.trial_type.sig_times(1,:) ==0)=0;
glm_output.pnc.sig_timePos = glm_output.pnc.sig_times & (glm_output.pnc.beta_weights>0);

% check that conflict modulation is positive
% glm_output.trial_type.sig_times(1,glm_output.trial_type.beta_weights(1,:)<0)=0;
% 
% %if noncacelled is higher, we don't care
% glm_output.nonCanc.sig_times(1,glm_output.nonCanc.beta_weights(1,:)<0) = 0;
% 

% 
% %now confirm that the SSD increases FR
% glm_output.ssd.sig_times(1,glm_output.ssd.beta_weights(1,:)<0)=0;
% 
% %only care about seen if SSD is significant
% glm_output.seen.sig_times(1,glm_output.ssd.sig_times(1,:)==0) = 0;
% 
% %only care about seen if seen is higher
% glm_output.seen.sig_times(1,glm_output.seen.beta_weights(1,:)<0) = 0;


signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 100ms of selectivity following SSD

% Trial-type
[canc_sig_start, canc_sig_len, ~] = ZeroOnesCount(glm_output.trial_type.sig_times(1,:)); % choice direction

% Stop-signal delay
[ssd_sig_start, ssd_sig_len, ~] = ZeroOnesCount(glm_output.ssd.sig_timePos(1,:) ); % choice direction

% inhibition function
[pnc_sig_start, pnc_sig_len, ~] = ZeroOnesCount(glm_output.pnc.sig_timePos(1,:) ); % choice direction

% nonCanc
[~, noncanc_sig_len, ~] = ZeroOnesCount(glm_output.nonCanc.sig_times(1,:)); % choice direction

% Seen
[~, seen_sig_len, ~] = ZeroOnesCount(glm_output.seen.sig_times(1,:)); % choice direction

% SSD perm
[~, ssdP_sig_len, ~] = ZeroOnesCount(glm_output.ssd.permSig(1,:)); % choice direction


encoding_flag = []; encoding_beta = [];


encoding_flag(1,1) = any(canc_sig_len >= signal_detect_wins);
encoding_flag(1,2) = any(ssd_sig_len >= signal_detect_wins);
encoding_flag(1,3) = any(pnc_sig_len >= signal_detect_wins);
encoding_flag(1,4) = any(noncanc_sig_len >= signal_detect_wins);
encoding_flag(1,5) = any(seen_sig_len >= signal_detect_wins);
encoding_flag(1,6) = any(ssdP_sig_len >= signal_detect_wins);
sig_bins = cell(1,3);
if encoding_flag(1,1)
    cancInd = find(canc_sig_len>=5,1);
    sig_bins{1,1} = canc_sig_start(cancInd):canc_sig_start(cancInd)+canc_sig_len(cancInd)-1;
    FRs = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
        [-60+(10*sig_bins{1}(1)) -60 + 100 + (10*sig_bins{1}(end))],canc_sig_len(cancInd)*10+90,window_shift);
    frs.canc{1} = mean(FRs(reg_tbl_trialtype.trial_number(reg_tbl_trialtype.trial_type==1)));
    frs.canc{2}= mean(FRs(reg_tbl_trialtype.trial_number(reg_tbl_trialtype.trial_type==2)));
end
pvals = struct();
if encoding_flag(1,2) 
    ssdInd = find(ssd_sig_len>=5,1);
    sig_bins{1,2} = ssd_sig_start(ssdInd):ssd_sig_start(ssdInd)+ssd_sig_len(ssdInd)-1;
    FRs = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
        [-60+(10*sig_bins{2}(1)) -60 + 100 + (10*sig_bins{2}(end))],ssd_sig_len(ssdInd)*10+90,window_shift);
    ssds = unique(reg_tbl_stopping.ssd);
    for ii = 1:length(ssds)
        latMatched = mean(FRs(reg_tbl_trialtype.trial_number(...
            reg_tbl_trialtype.ssd == ssds(ii) & reg_tbl_trialtype.trial_type ==1)));
        frs.ssd{ii} = mean(FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))));
        tempfrs{ii} = FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))) - latMatched;
        tempssd{ii}  = ones(length(tempfrs{ii}),1)*inh_pnc(2,ssds(ii));
        tempPNC{ii} = ones(length(tempfrs{ii}),1)*inh_pnc(1,ssds(ii));
        frs.ssdDiff{ii} = frs.ssd{ii}-latMatched;
    end
    ssdfitssd = fitlm(vertcat(tempssd{:}), vertcat(tempfrs{:}));
    ssdfitssdp = ssdfitssd.anova.pValue;
    ssdfitpnc = fitlm(vertcat(tempPNC{:}), vertcat(tempfrs{:}));
    ssdfitpncp = ssdfitpnc.anova.pValue;
    pvals.SSD = table(ssdfitssdp,ssdfitpncp);
end
if encoding_flag(1,3)
    pncInd = find(pnc_sig_len>=5,1);
    sig_bins{1,3} = pnc_sig_start(pncInd):pnc_sig_start(pncInd)+pnc_sig_len(pncInd)-1;

    FRs = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
        [-60+(10*sig_bins{3}(1)) -60 + 100 + (10*sig_bins{3}(end))],pnc_sig_len(pncInd)*10+90,window_shift);
    ssds = unique(reg_tbl_stopping.ssd);
    for ii = 1:length(ssds)
         latMatched = mean(FRs(reg_tbl_trialtype.trial_number(...
            reg_tbl_trialtype.ssd == ssds(ii) & reg_tbl_trialtype.trial_type ==1)));
      
        frs.pnc{ii} = mean(FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))));
        frs.pncDiff{ii} = frs.pnc{ii} - latMatched;
         tempfrs{ii} = FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))) - latMatched;
        tempssd{ii}  = ones(length(tempfrs{ii}),1)*inh_pnc(2,ssds(ii));
        tempPNC{ii} = ones(length(tempfrs{ii}),1)*inh_pnc(1,ssds(ii));
    end
    pncfitssd = fitlm(vertcat(tempssd{:}), vertcat(tempfrs{:}));
    pncfitssdp = pncfitssd.anova.pValue;
    pncfitpnc = fitlm(vertcat(tempPNC{:}), vertcat(tempfrs{:}));
    pncfitpncp = pncfitpnc.anova.pValue;
    pvals.PNC = table(pncfitssdp,pncfitpncp);
end
% Get average beta-weights
encoding_beta(1,1) = nanmean(glm_output.trial_type.beta_weights(1,sig_bins{1}));
encoding_beta(1,2) = nanmean(glm_output.ssd.beta_weights(1,sig_bins{2}));
encoding_beta(1,3) = nanmean(glm_output.pnc.beta_weights(1,sig_bins{3}));
encoding_beta(1,4) = nanmean(glm_output.nonCanc.beta_weights(1,...
    glm_output.nonCanc.sig_times(1,:)==1));
encoding_beta(1,5) = nanmean(glm_output.seen.beta_weights(1,...
    glm_output.seen.sig_times(1,:)==1));



%% Workpad: figure checks
% % 
% fr_times = -0:10:450;
% window_time = -1000:2000;
% figure;
% subplot(4,2,[1 3]); hold on
% plot(-1000:2000,nanmean(sdfs(reg_tbl.trial_type==2,:)),'r')
% plot(window_time,nanmean(sdfs(reg_tbl.trial_type==1,:)),'k')
% xlim([-250 1000])
% % 
% subplot(4,2,5)
% imagesc('XData',fr_times,'YData',ones(1,length(window_time)),'CData',glm_output.trial_type.sig_times(1,:))
% xlim([-250 1000])
% 
% subplot(4,2,7)
% plot(fr_times, glm_output.trial_type.beta_weights(1,:))
% xlim([-250 1000])
% 
% 
% subplot(4,2,[2 4]); hold on
% plot(window_time,nanmean(sdfs(reg_tbl.ssd== 1,:)),'color',[1 0 0 0.33])
% plot(window_time,nanmean(sdfs(reg_tbl.ssd==2,:)),'color',[1 0 0 0.66])
% plot(window_time,nanmean(sdfs(reg_tbl.ssd==3,:)),'color',[1 0 0 0.99])
% xlim([-250 1000])
% 
% subplot(4,2,6)
% imagesc('XData',fr_times,'YData',ones(1,length(fr_times)),'CData',glm_output.ssd.sig_times(1,:))
% xlim([-250 1000])
% 
% subplot(4,2,8)
% plot(fr_times, glm_output.ssd.beta_weights(1,:))
% xlim([-250 1000])
% 
% 
% 
% 
% 
% 
% 
% 
% 
