
%% Neuron GLM 
% Approx 2hr run time; output saved 2023-09-03, 20h57 (stopping_glm_out)
files = dir('E:\conflict\results');
fileNames = {files(:).name};
fileCheck = cell2mat(cellfun(@(x) contains(x,[region 'SDFDiff_workspaceA102dirPNC']),fileNames,'uni',0));
if any(fileCheck)
    load(['E:\conflict\results\' fileNames{fileCheck}])
else
    rng("default")
parfor neuron_i = 1:size(mcc_map_info,1)
    neuralFilename = mcc_map_info.session{neuron_i};%... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
    warning off
    [glm_out_mcc{neuron_i,1}, encoding_flag_mcc(neuron_i,:), ...
        inh_pNC{neuron_i,1}] =...
        stopping_diffSDF_fr_aligned(neuron_i, mcc_map_info, behavior(behaviorIdx,:));
end
% 
% parfor neuron_i = 1:size(mcc_map_info,1)
%     neuralFilename = mcc_map_info.session{neuron_i};
%     behFilename = data_findBehFile(neuralFilename);
%     behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
%     neuronLabel = mcc_map_info.unit{neuron_i};
%     fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
%     warning off
%     [glm_out_mccB{neuron_i,1}, encoding_flag_mccB(neuron_i,:), encoding_beta_mccB(neuron_i,:)] =...
%         stopping_glm_func_fr_aligned_shuffle(neuron_i, mcc_map_info, behavior(behaviorIdx,:));
% end
% parfor neuron_i = 1:size(mcc_map_info,1)
%     neuralFilename = mcc_map_info.session{neuron_i};
%      behFilename = data_findBehFile(neuralFilename);
%     behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
%     neuronLabel = mcc_map_info.unit{neuron_i};
%     fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
%     warning off
%     [glm_out_mccB{neuron_i,1}, encoding_flag_mccB(neuron_i,:), encoding_beta_mccB(neuron_i,:)] =...
%         stopping_glm_func_fr_bline(neuron_i, mcc_map_info, dataFiles_beh, dirs,behavior(behaviorIdx,:));
% end

%% Extract stop-related spike density functions
parfor neuron_i = 1:size(mcc_map_info,1)
    
    neuralFilename = mcc_map_info.session{neuron_i};
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    %load in stationarity data
    statData = load(['E:\conflict\results\DATA_SPK_STATIONARITY\Ben\stationarityData_' mcc_map_info.probeName{neuron_i}])
    ndata = load(['E:\conflict\results\probeTrls\neuronData_' mcc_map_info.probeName{neuron_i}])
    tempInd = find(contains(ndata.neuronData.DSPname,neuronLabel));
    stationInd = find(statData.stationaritySelection.included_trials.neuronIdx == tempInd);
    goodTrls = statData.stationaritySelection.included_trials.trials{stationInd};
    % Load in pre-processed spike data
    data_in = load(fullfile('E:\conflict\data\','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    event_alignment = 'target';
    baseline_win = [-600:-100];
    
    z_sdf = zscore_sdf(data_in.SDF,baseline_win,event_alignment);
    allTrls = true(size(z_sdf,1),1);   allTrls(goodTrls) = false;
    z_sdf(allTrls,:) = nan(sum(allTrls==true),3001);
    SSDs{neuron_i} = [behavior(behaviorIdx,:).stopSignalBeh.inh_SSD,...
        behavior(behaviorIdx,:).stopSignalBeh.inh_pnc'];
    ssrt_delays = round(behavior(behaviorIdx,:).stopSignalBeh.inh_SSD + ...
        behavior(behaviorIdx,:).stopSignalBeh.ssrt.integrationWeighted);
    sdf_canceled_ssdx = []; sdf_nostop_ssdx = []; sdf_noncanc_ssdx = [];
    sdf_seen_ssdx = [];     sdf_unseen_ssdx = [];

    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx,:).stopSignalBeh.inh_SSD)
        trl_canceled = []; trl_canceled = behavior(behaviorIdx,:).ttm.C.C{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx,:).ttm.C.GO{ssd_i};
        trl_noncanc = []; trl_noncanc = behavior(behaviorIdx,:).ttm.NC.NC{ssd_i};
        trl_seen = []; trl_seen = behavior(behaviorIdx,:).ttseen.NC.S{ssd_i};
        trl_unseen = []; trl_unseen = behavior(behaviorIdx,:).ttseen.NC.U{ssd_i};
        
        trlLim = 7;
        if length(trl_canceled) < trlLim || length(trl_nostop) < trlLim ||...
                length(trl_noncanc) < 2
            sdf_canceled_ssdx(ssd_i,:) = nan(1,length(-1000:1000));
            sdf_nostop_ssdx(ssd_i,:) = nan(1,length(-1000:1000));
            sdf_noncanc_ssdx(ssd_i,:) = nan(1,length(-1000:1000));
            sdf_seen_ssdx(ssd_i,:) = nan(1,length(-1000:1000));
            sdf_unseen_ssdx(ssd_i,:) = nan(1,length(-1000:1000));

        else
            %because we are using target aligned sdfs from -1000:2000, we
            %can realign and get -1000:1000 by shifting the start by the
            %SSD + SSRT and then capture the next 2000 ms. now 1000 is
            %aligned to SSRT
            sdf_canceled_ssdx(ssd_i,:) = nanmean(z_sdf(trl_canceled,ssrt_delays(ssd_i):...
                ssrt_delays(ssd_i) + 2000));
            sdf_nostop_ssdx(ssd_i,:) = nanmean(z_sdf(trl_nostop,ssrt_delays(ssd_i):...
                ssrt_delays(ssd_i) + 2000));      
            sdf_noncanc_ssdx(ssd_i,:) = nanmean(z_sdf(trl_noncanc,ssrt_delays(ssd_i):...
                ssrt_delays(ssd_i) + 2000));      
            sdf_seen_ssdx(ssd_i,:) = nanmean(z_sdf(trl_seen,ssrt_delays(ssd_i):...
                ssrt_delays(ssd_i) + 2000));      
            sdf_unseen_ssdx(ssd_i,:) = nanmean(z_sdf(trl_unseen,ssrt_delays(ssd_i):...
                ssrt_delays(ssd_i) + 2000),1);      
        end
    end
    
    diff_cancelled{neuron_i,:} = sdf_canceled_ssdx - sdf_nostop_ssdx;
    sdf_canceled{neuron_i,:} = sdf_canceled_ssdx;
    sdf_nostop{neuron_i,:} = sdf_nostop_ssdx;
    sdf_noncanc{neuron_i,:} = sdf_noncanc_ssdx;
    sdf_seen{neuron_i,:} = sdf_seen_ssdx;
    sdf_unseen{neuron_i,:} = sdf_unseen_ssdx;
    seenCounts = cellfun(@length, behavior(behaviorIdx,:).ttseen.NC.S);
    unseenCounts = cellfun(@length, behavior(behaviorIdx,:).ttseen.NC.U);
    sdf_seenAll(neuron_i,:) = nansum(sdf_seen_ssdx .* seenCounts')./sum(seenCounts);
    sdf_unseenAll(neuron_i,:) = nansum(sdf_unseen_ssdx .* unseenCounts')./sum(unseenCounts);
     NScounts = cellfun(@length,behavior(behaviorIdx,:).ttm.C.GO);
    cancCounts = cellfun(@length,behavior(behaviorIdx,:).ttm.C.C);
    sdf_canceledAll(neuron_i,:) = nansum(sdf_canceled_ssdx .* cancCounts')./sum(cancCounts);
    sdf_nostopAll(neuron_i,:) = nansum(sdf_nostop_ssdx .* NScounts')./sum(NScounts);
    
    
    sdf_highRwd(neuron_i,:) = nanmean(z_sdf(behavior(behaviorIdx,:).ttx.canceled.all.hi,:));
    sdf_lowRwd(neuron_i,:) = nanmean(z_sdf(behavior(behaviorIdx,:).ttx.canceled.all.lo,:));
end
    save(['E:\conflict\results\' region 'SDFDiff_workspaceA102dirPNC.mat'])

end
%% Export: export table for future analyses
mcc_analysis_table = table();
mcc_analysis_table.index = [1:size(mcc_map_info,1)]';
mcc_analysis_table.session = mcc_map_info.session;
mcc_analysis_table.unit = mcc_map_info.unit;

mcc_analysis_table.sdf_canceledDiff = diff_cancelled;
mcc_analysis_table.sdf_canceled = sdf_canceled;
mcc_analysis_table.sdf_nostop = sdf_nostop;
mcc_analysis_table.sdf_noncanc = sdf_noncanc;
mcc_analysis_table.sdf_seen = sdf_seen;
mcc_analysis_table.sdf_unseen = sdf_unseen;
mcc_analysis_table.sdf_seenAll = sdf_seenAll;
mcc_analysis_table.sdf_unseenAll = sdf_unseenAll;
mcc_analysis_table.sdf_canceledAll = sdf_canceledAll;
mcc_analysis_table.sdf_nostopAll = sdf_nostopAll;

mcc_analysis_table.sdf_highRwd = sdf_highRwd;
mcc_analysis_table.sdf_lowRwd = sdf_lowRwd;

mcc_analysis_table.SSDs = SSDs(:,1);

mcc_analysis_table.glm_trial = encoding_flag_mcc(:,1);
mcc_analysis_table.glm_trialPos = encoding_flag_mcc(:,1) &encoding_beta_mcc(:,1)>0;
mcc_analysis_table.glm_trialNeg = encoding_flag_mcc(:,1) &encoding_beta_mcc(:,1)<0;

mcc_analysis_table.glm_ssd_canc = encoding_flag_mcc(:,2) ;
mcc_analysis_table.glm_ssd_cancPos = encoding_flag_mcc(:,2) & encoding_beta_mcc(:,2)>0;
mcc_analysis_table.glm_ssd_cancNeg = encoding_flag_mcc(:,2) & encoding_beta_mcc(:,2)<0;

mcc_analysis_table.glm_pnc_canc = encoding_flag_mcc(:,3) ;
mcc_analysis_table.glm_pnc_cancPos = encoding_flag_mcc(:,3) & encoding_beta_mcc(:,3)>0;
mcc_analysis_table.glm_pnc_cancNeg = encoding_flag_mcc(:,3) & encoding_beta_mcc(:,3)<0;

mcc_analysis_table.glm_ssd_cancP = encoding_flag_mcc(:,5) ;
mcc_analysis_table.glm_ssd_cancPosP = encoding_flag_mcc(:,5) & encoding_beta_mcc(:,2)>0;
mcc_analysis_table.glm_ssd_cancNegP = encoding_flag_mcc(:,5) & encoding_beta_mcc(:,2)<0;

mcc_analysis_table.glm_non_canc = encoding_flag_mcc(:,3);
mcc_analysis_table.glm_seen = encoding_flag_mcc(:,4);

%% Index GLM +ve neurons
trial_type_neurons = [];
trial_type_neurons = find(mcc_analysis_table.glm_ssd_canc == 1);

%% Find modulation direction
glm_times = [0:10:300];
glm_index_ref = [1:31];
modulation_window = [0:600];
canc_nostop_fr_diff = zeros(size(trial_type_neurons,1),1);
for neuron_i = 1:size(trial_type_neurons,1)
    neuron_j = trial_type_neurons(neuron_i);
    
    glm_window = [];
    glm_window =...
        glm_times(find(glm_out_mcc{neuron_j}.trial_type.sig_times(1,glm_index_ref) == 1,1,'first')):...
        glm_times(find(glm_out_mcc{neuron_j}.trial_type.sig_times(1,glm_index_ref) == 1,1,'last'));
    
    canc_fr = nanmean(nanmean(mcc_analysis_table.sdf_canceled{neuron_j}(:,glm_window+1000)));
    nostop_fr = nanmean(nanmean(mcc_analysis_table.sdf_nostop{neuron_j}(:,glm_window+1000)));
    canc_nostop_fr_diff(neuron_i) = (canc_fr-nostop_fr)/abs(nostop_fr);
    
end

% Plot histogram of % diff
figure;
histogram(canc_nostop_fr_diff,0:0.025:5,'LineStyle','None'); %xlim([0 2]); 
vline(1)
xlabel('Canceled firing % relative to no-stop'); ylabel('Count')

%% Place neurons into structure divided by type
neuron_index.trial_type.nostop = trial_type_neurons(canc_nostop_fr_diff  <= -0.1);
neuron_index.trial_type.canceled = trial_type_neurons(canc_nostop_fr_diff  >= 0.1);
neuron_index.ssd.facAscend = find(mcc_analysis_table.glm_ssd_cancPos...
    & mcc_analysis_table.glm_trialPos);
neuron_index.ssd.facDescend = find(mcc_analysis_table.glm_ssd_cancNeg...
    & mcc_analysis_table.glm_trialPos);
neuron_index.ssd.supAscend = find(mcc_analysis_table.glm_ssd_cancPos...
    & mcc_analysis_table.glm_trialNeg);
neuron_index.ssd.supDescend = find(mcc_analysis_table.glm_ssd_cancNeg...
    & mcc_analysis_table.glm_trialNeg);

neuron_index.pnc.facAscend = find(mcc_analysis_table.glm_pnc_cancPos...
    & mcc_analysis_table.glm_trialPos);
neuron_index.pnc.facDescend = find(mcc_analysis_table.glm_pnc_cancNeg...
    & mcc_analysis_table.glm_trialPos);
neuron_index.pnc.supAscend = find(mcc_analysis_table.glm_pnc_cancPos...
    & mcc_analysis_table.glm_trialNeg);
neuron_index.pnc.supDescend = find(mcc_analysis_table.glm_pnc_cancNeg...
    & mcc_analysis_table.glm_trialNeg);

neuron_index.trial_type.facAscendP = find(mcc_analysis_table.glm_ssd_cancPosP...
    & mcc_analysis_table.glm_trialPos);
neuron_index.trial_type.facDescendP = find(mcc_analysis_table.glm_ssd_cancNegP...
    & mcc_analysis_table.glm_trialPos);
neuron_index.trial_type.supAscendP = find(mcc_analysis_table.glm_ssd_cancPosP...
    & mcc_analysis_table.glm_trialNeg);
neuron_index.trial_type.supDescendP = find(mcc_analysis_table.glm_ssd_cancNegP...
    & mcc_analysis_table.glm_trialNeg);
% neuron_index.trial_type.nostop = trial_type_neurons(canc_nostop_fr_diff  <= 0.9);
% neuron_index.trial_type.canceled = trial_type_neurons(canc_nostop_fr_diff  >= 1.1);
alpha = .02;
FacSup = 'Fac';
[output,permEnc_flag,permEnc_beta] = permAssessment(glm_out_mcc,FacSup,alpha); %alpha is two sided, so double desired alpha for 1-sided
neuron_index.trial_type.permAscFac = find(permEnc_flag(:,2));
neuron_index.trial_type.permTrialTypeFac = find(permEnc_flag(:,1));
neuron_index.trial_type.permDescFac = find(permEnc_flag(:,3));
FacSup = 'Sup';
[output,permEnc_flag,permEnc_beta] = permAssessment(glm_out_mcc,FacSup,alpha); %alpha is two sided, so double desired alpha for 1-sided
neuron_index.trial_type.permAscSup = find(permEnc_flag(:,2));
neuron_index.trial_type.permTrialTypeSup = find(permEnc_flag(:,1));
neuron_index.trial_type.permDescSup = find(permEnc_flag(:,3));


mcc_analysis_table.canc_type = ismember(1:size(mcc_analysis_table,1),neuron_index.trial_type.canceled)';
mcc_analysis_table.nostop_type = ismember(1:size(mcc_analysis_table,1),neuron_index.trial_type.nostop)';



% %% Figure: produce a venn diagram of the counts (canceled)
% mcc_analysis_table_canceled = [];
% mcc_analysis_table_canceled = mcc_analysis_table(neuron_index.trial_type.canceled,:);
% 
% glm_counts_canc = [];
% glm_counts_canc = ...
%     [... % One factor
%     sum(mcc_analysis_table_canceled.glm_trial == 1 & mcc_analysis_table_canceled.glm_ssd_canc == 0 & mcc_analysis_table_canceled.glm_non_canc == 0);...
%     sum(mcc_analysis_table_canceled.glm_trial == 0 & mcc_analysis_table_canceled.glm_ssd_canc == 1 & mcc_analysis_table_canceled.glm_non_canc == 0);...
%     sum(mcc_analysis_table_canceled.glm_trial == 0 & mcc_analysis_table_canceled.glm_ssd_canc == 0 & mcc_analysis_table_canceled.glm_non_canc == 1);...
%     % Two factors
%     sum(mcc_analysis_table_canceled.glm_trial == 1 & mcc_analysis_table_canceled.glm_ssd_canc == 1 & mcc_analysis_table_canceled.glm_non_canc == 0);...
%     sum(mcc_analysis_table_canceled.glm_trial == 1 & mcc_analysis_table_canceled.glm_ssd_canc == 0 & mcc_analysis_table_canceled.glm_non_canc == 1);...
%     sum(mcc_analysis_table_canceled.glm_trial == 0 & mcc_analysis_table_canceled.glm_ssd_canc == 1 & mcc_analysis_table_canceled.glm_non_canc == 1);...
%     % All factors
%     sum(mcc_analysis_table_canceled.glm_trial == 1 & mcc_analysis_table_canceled.glm_ssd_canc == 1 & mcc_analysis_table_canceled.glm_non_canc == 1)];
% 
% 
% mysets = ["Action" "SSD" "Value"];
% mylabels = glm_counts_canc;
% venn(3,'sets',mysets,'edgeC',[0 0 0],'colors',autumn(3),'labels',mylabels,'edgeW',2);
% sum(glm_counts_canc);
% 
% %% Figure: produce a venn diagram of the counts (nostop)
% mcc_analysis_table_nostop = [];
% mcc_analysis_table_nostop = mcc_analysis_table(neuron_index.trial_type.nostop,:);
% 
% glm_counts_nostop = [];
% glm_counts_nostop = ...
%     [... % One factor
%     sum(mcc_analysis_table_nostop.glm_trial == 1 & mcc_analysis_table_nostop.glm_ssd_canc == 0 & mcc_analysis_table_nostop.glm_value_canc == 0);...
%     sum(mcc_analysis_table_nostop.glm_trial == 0 & mcc_analysis_table_nostop.glm_ssd_canc == 1 & mcc_analysis_table_nostop.glm_value_canc == 0);...
%     sum(mcc_analysis_table_nostop.glm_trial == 0 & mcc_analysis_table_nostop.glm_ssd_canc == 0 & mcc_analysis_table_nostop.glm_value_canc == 1);...
%     % Two factors
%     sum(mcc_analysis_table_nostop.glm_trial == 1 & mcc_analysis_table_nostop.glm_ssd_canc == 1 & mcc_analysis_table_nostop.glm_value_canc == 0);...
%     sum(mcc_analysis_table_nostop.glm_trial == 1 & mcc_analysis_table_nostop.glm_ssd_canc == 0 & mcc_analysis_table_nostop.glm_value_canc == 1);...
%     sum(mcc_analysis_table_nostop.glm_trial == 0 & mcc_analysis_table_nostop.glm_ssd_canc == 1 & mcc_analysis_table_nostop.glm_value_canc == 1);...
%     % All factors
%     sum(mcc_analysis_table_nostop.glm_trial == 1 & mcc_analysis_table_nostop.glm_ssd_canc == 1 & mcc_analysis_table_nostop.glm_value_canc == 1)];
% 
% mysets = ["Action" "SSD" "Value"];
% mylabels = glm_counts_nostop;
% venn(3,'sets',mysets,'edgeC',[0 0 0],'colors',summer(3),'labels',mylabels,'edgeW',2);
% sum(glm_counts_nostop);


