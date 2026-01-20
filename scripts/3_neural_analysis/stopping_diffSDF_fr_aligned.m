function [glm_output, encoding_flag, inh_pnc,errTable] = ...
    stopping_diffSDF_fr_aligned(neuron_i, mcc_map_info, behData,dirs)
rng(99999999)
randi(50)
neuralFilename = mcc_map_info.session{neuron_i};
neuronLabel = mcc_map_info.unit{neuron_i};
% fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
 

% Load in pre-processed spike data 
% import_data = struct(); import_data = load_behFile(dirs,behFilename);
regAlpha = .05;
doPerm = 0;
contFiller = 20;
load([ dirs.results '\DATA_SPK_STATIONARITY\Ben\stationarityData_' mcc_map_info.probeName{neuron_i}])
load([ dirs.results '\probeTrls\neuronData_' mcc_map_info.probeName{neuron_i}])
tempInd = find(contains(neuronData.DSPname,neuronLabel));
stationInd = stationaritySelection.included_trials.neuronIdx == tempInd;
rightInds = getDir(behData.ttx);

%% Extract: get relevant data for GLM table

reg_tbl = table;

reg_tbl.trial_type = zeros(length(behData.trialEventTimes{1}.isNoStop),1);
reg_tbl.trial_type(behData.trialEventTimes{1}.isNoStop==1) = 1;
reg_tbl.trial_type(behData.trialEventTimes{1}.isCanceled==1) = 2;
reg_tbl.trial_type(behData.trialEventTimes{1}.isNonCanceled==1) = 3;
reg_tbl.seen(~isnan(behData.trialEventTimes{1}.stopSignal)) = 1;
reg_tbl.rewarded(~isnan(behData.trialEventTimes{1}.reward) & ...
    reg_tbl.trial_type ~=3) = 1;
reg_tbl.ssd = behData.trialEventTimes{1}.ssdInd+1;
reg_tbl.SS = round(behData.trialEventTimes{1}.stopSignal - behData.trialEventTimes{1}.target);
reg_tbl.RT = round(behData.trialEventTimes{1}.saccade - behData.trialEventTimes{1}.target);
reg_tbl.errRT = round(behData.trialEventTimes{1}.saccade - behData.trialEventTimes{1}.stopSignal);
reg_tbl.seen(reg_tbl.errRT<1) = 0; %change all unseen stop signals to seen==0
reg_tbl.trial_number = [1:height(reg_tbl)]';
reg_tbl.isR = false(height(reg_tbl),1);
reg_tbl.isR(rightInds) = true;
reg_tbl.isLow = behData.trialEventTimes{1}.isLowRwd;
reg_tbl = expectationRew(reg_tbl);


goodtrls = false(height(reg_tbl),1);
goodtrls(stationaritySelection.included_trials.trials{stationInd}) = true;
%map ssd to surprise
surprise = zeros(1,20);
for jj = 1:length(behData.stopSignalBeh.inh_SSD)
for ii = 1:length(behData.Hazard{1})
 if behData.stopSignalBeh.inh_SSD(jj) > (40*ii) & behData.stopSignalBeh.inh_SSD(jj) <= (40 * (ii+1 ))
surprise(jj) = -log2(behData.Hazard{1}(ii));
 end
end

end
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
    reg_tbl.surprise(reg_tbl.ssd==ii) = surprise(ii);
    reg_tbl.ssdVal(reg_tbl.ssd==ii) = behData.stopSignalBeh.inh_SSD(ii);
    inh_pnc(1,ii) = behData.stopSignalBeh.inh_pnc(ii);
    inh_pnc(2,ii) = behData.stopSignalBeh.inh_SSD(ii)';
    inh_pnc(3,ii) = surprise(ii);
end
%intialize fr struct
frs.canc = {};
frs.ssd = {};
frs.pnc = {};
if height(reg_tbl(goodtrls,:))==0 || sum(reg_tbl.trial_type(goodtrls)==2)<30
    glm_output.trial_typeSDF.CIs = [];
    glm_output.trial_typeSDF.diff = [];
    glm_output.trial_typeSDF.table = [];
    glm_output.trial_canNCSDF.CIs = [];
    glm_output.trial_canNCSDF.diff = [];
    glm_output.trial_canNCSDF.table = [];
    glm_output.trial_errorGoSDF.CIs = [];
    glm_output.trial_errorGoSDF.diff = [];
    glm_output.trial_errorGoSDF.table = [];
    glm_output.trial_errorGoSDF.cancSDF = [];
    glm_output.trial_errorGoSDF.noCancSDF = [];
    glm_output.trial_errorGoSDFCorr.CIs = [];
    glm_output.trial_errorGoSDFCorr.diff = [];
    glm_output.trial_errorGoSDFCorr.table = [];
    glm_output.trial_errorGoSDFCorr.cancSDF = [];
    glm_output.trial_errorGoSDFCorr.noCancSDF = [];
    glm_output.trial_errorGoSDFCorr.AUC = [];
    glm_output.trial_errorGoSDFCorr.AUC2sd = [];
    glm_output.trial_errorGoSDFCorr.AUCNarr = [];
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
    encoding_flag = zeros(1,7); 
    sig_bins = cell(1,3);
    pvals = struct();
    remove = 1;
    return
end
%% Setup spike data into GLM
spk_data_in = load(fullfile(dirs.data,'Spikes',...
    [neuralFilename '_Spikes_' neuronLabel '.mat']));
event_alignment = 'ssrt';
% allSDFSS = getIndividSDF_MCC(spk_data_in.Spikes,event_alignment);
window_shift = 10;
fr_window = 100;
baseline_win = [-600:-500];

%% Run GLM: trial type
glm_output.trial_type.sig_times = []; 
glm_output.trial_type.beta_weights = []; 
glm_output.trial_type.permPos = []; 
glm_output.trial_type.pvals = []; 
u_t_mdl = [];
reg_tbl.is_canc = reg_tbl.trial_type==2;
reg_tbl.trl_numZ = zscore(reg_tbl.trial_number);
%match trials to calculate noStop firing rates
SSRT = behData.stopSignalBeh.ssrt{3}.integrationWeighted;
ssrt_delays = round(behData.stopSignalBeh.inh_SSD + ...
    SSRT);

trlInds = reg_tbl.trial_number;
goodTrlNums = find(goodtrls);
[matchedInds,matchedBySSD] = trlMatcherCanc(behData,goodTrlNums,reg_tbl);
[~,event_SDFs,eventRasters] = calc_FRs_matched(spk_data_in.Spikes,baseline_win,event_alignment,...
    ssrt_delays,matchedBySSD,[-50, 350],fr_window,window_shift);
% reg_tbl.win_fr = event_frs_matched;
reg_tbl.sdfs = event_SDFs;
reg_tbl.matched_rasters = eventRasters;
%only used for the SSD scaling on stop trials
reg_tbl.rasters = spk_data_in.Spikes.ssrt;
reg_tbl.rastersSac = spk_data_in.Spikes.saccade;
%assign the matched SSD
for ssd_i = 1:length(matchedBySSD)
    reg_tbl.ssd(matchedBySSD{ssd_i,1}) = ssd_i;
end
matchedInds = trlInds(ismember(trlInds,matchedInds));

reg_tbl_trialtype = reg_tbl(matchedInds,:);
cancelledSDF = mean(reg_tbl_trialtype.sdfs(reg_tbl_trialtype.trial_type==2,:));
noStopSDF = mean(reg_tbl_trialtype.sdfs(reg_tbl_trialtype.trial_type==1,:));
[cancGoCI, cancGoDiff, ~,actualData] = ...
                calculatePermutationContinuity(reg_tbl_trialtype.sdfs,...
                reg_tbl_trialtype.trial_type==2,reg_tbl_trialtype.trial_type==1,contFiller);
glm_output.trial_typeSDF.CIs = cancGoCI;
glm_output.trial_typeSDF.diff = cancGoDiff;
glm_output.trial_typeSDF.table = actualData;
glm_output.trial_typeSDF.cancSDF = cancelledSDF;
glm_output.trial_typeSDF.nostopSDF = noStopSDF;
%% Cancelled vs error

[matchedInds,matchedBySSDNC] = trlMatcherCancNC(behData,goodTrlNums,reg_tbl);
[~,event_SDFs] = calc_FRs_matched(spk_data_in.Spikes,baseline_win,event_alignment,...
    ssrt_delays,matchedBySSDNC,[-50, 350],fr_window,window_shift);
% reg_tbl.win_fr = event_frs_matched;
matchedInds = trlInds(ismember(trlInds,matchedInds));

reg_tbl.sdfs(matchedInds,:) = event_SDFs(matchedInds,:);

reg_tbl_canc = reg_tbl(matchedInds,:);
cancelledSDF = mean(reg_tbl_canc.sdfs(reg_tbl_canc.trial_type==2,:));
noncancelledSDF = mean(reg_tbl_canc.sdfs(reg_tbl_canc.trial_type==3,:));

[cancNCCI, cancNCDiff, ~,actualDataNC] = ...
                calculatePermutationContinuity(reg_tbl_canc.sdfs,...
                reg_tbl_canc.trial_type==2,reg_tbl_canc.trial_type==3,contFiller);
cancelledSDFNC = mean(reg_tbl_canc.sdfs(reg_tbl_canc.trial_type==2,:));
nonCancSDF = mean(reg_tbl_canc.sdfs(reg_tbl_canc.trial_type==3,:));

glm_output.trial_canNCSDF.CIs = cancNCCI;
glm_output.trial_canNCSDF.diff = cancNCDiff;
glm_output.trial_canNCSDF.table = actualDataNC;
glm_output.trial_canNCSDF.cancSDF = cancelledSDF;
glm_output.trial_canNCSDF.noCancSDF = noncancelledSDF;

%% error vs nostop
[matchedInds,matchedBySSDSeen] = trlMatcherNC_proximity(behData,goodTrlNums,reg_tbl);

[~,event_SDFs,eventRasters] = calc_FRs_matched_saccade(spk_data_in.Spikes,baseline_win,'saccade',...
    matchedBySSDSeen,[-50, 350],fr_window,window_shift);
% reg_tbl.win_fr = event_frs_matched;
matchedInds = trlInds(ismember(trlInds,matchedInds));

reg_tbl.sdfs3k(matchedInds,:) = event_SDFs(matchedInds,:);
reg_tbl.matched_rasters3k(matchedInds,:) = eventRasters(matchedInds,1:3001);

reg_tbl_error = reg_tbl(matchedInds,:);
rng(999)
[errorCI, errorDiff, ~,actualDataError] = ...
                calculatePermutationContinuity(reg_tbl_error.sdfs3k,...
                reg_tbl_error.trial_type==1,reg_tbl_error.trial_type==3,contFiller);
noStopSDFErr = mean(reg_tbl_error.sdfs3k(reg_tbl_error.trial_type==1,:));
nonCancSDFErr = mean(reg_tbl_error.sdfs3k(reg_tbl_error.trial_type==3,:));


bfr = mean(sum(eventRasters(matchedInds,1:501),2)./.5);
glm_output.trial_errorGoSDF.CIs = errorCI;
glm_output.trial_errorGoSDF.diff = errorDiff;
glm_output.trial_errorGoSDF.table = actualDataError;
glm_output.trial_errorGoSDF.nostopSDF = noStopSDFErr;
glm_output.trial_errorGoSDF.noCancSDF = nonCancSDFErr;
glm_output.trial_errorGoSDF.baselineFR = bfr;

%% error by value
% figure
% tiledlayout
% nexttile
% hold on
% plot(-1000:2000,mean(event_SDFs(matchedBySSDSeen{1},:))*1000)
% plot(-1000:2000,mean(event_SDFs(matchedBySSDSeen{2},:))*1000)
% xlim([-200 600])
% legend({'no stop','error'})
goodTrlErrors = find(goodtrls & (reg_tbl.lHi | reg_tbl.rHi));
[matchedIndsTemp,matchedBySSDSeenTemp] = trlMatcherError_noRT(behData,goodTrlErrors,reg_tbl);
[matchedIndsCorr,matchedBySSDSeenCorr] = trlMatcherCorr_noRT(behData,goodTrlErrors,reg_tbl);
[~,event_SDFs,eventRasterstemp] = calc_FRs_matched_saccade(spk_data_in.Spikes,baseline_win,'saccade',...
    matchedBySSDSeenTemp,[-50, 350],fr_window,window_shift);
[~,event_SDFsCorr,eventRasterstempCorr] = calc_FRs_matched_saccade(spk_data_in.Spikes,baseline_win,'saccade',...
    matchedBySSDSeenCorr,[-50, 350],fr_window,window_shift);
eventRasterstempCorr(matchedIndsTemp,:) = eventRasterstemp(matchedIndsTemp,:);
event_SDFsCorr(matchedIndsTemp,:) = event_SDFs(matchedIndsTemp,:);
reg_tbl.matched_rasters3k(sort([matchedIndsTemp;matchedIndsCorr]),:) =...
    eventRasterstempCorr(sort([matchedIndsTemp;matchedIndsCorr]),:);
reg_tbl.sdfs3k(sort([matchedIndsTemp;matchedIndsCorr]),:) =...
    event_SDFsCorr(sort([matchedIndsTemp;matchedIndsCorr]),:);

reg_tbl_hilo = reg_tbl(sort([matchedIndsTemp;matchedIndsCorr]),:);

[errorCI, errorDiff, ~,actualDataHiLo] = ...
    calculatePermutationContinuity(reg_tbl_hilo.sdfs3k,...
    reg_tbl_hilo.trial_type==3 & reg_tbl_hilo.isH,....
    reg_tbl_hilo.trial_type==3 & ~reg_tbl_hilo.isH,20);


ErrHi = mean(reg_tbl_hilo.sdfs3k(reg_tbl_hilo.trial_type==3 & reg_tbl_hilo.isH,:));
errLo = mean(reg_tbl_hilo.sdfs3k(reg_tbl_hilo.trial_type==3 & ~reg_tbl_hilo.isH,:));
corrHi = mean(reg_tbl_hilo.sdfs3k(reg_tbl_hilo.trial_type==1 & reg_tbl_hilo.isH,:));
corrLo = mean(reg_tbl_hilo.sdfs3k(reg_tbl_hilo.trial_type==1 & ~reg_tbl_hilo.isH,:));

frs = sum(eventRasters(:,1000:1200),2);
% nexttile
% hold on
% scatter(reg_tbl.RT(matchedBySSDSeen{1}),frs(matchedBySSDSeen{1}))
% scatter(reg_tbl.RT(matchedBySSDSeen{2}),frs(matchedBySSDSeen{2}))
[rall,p] = corrcoef(reg_tbl.RT(matchedInds), frs(matchedInds));
% title(['r = ' num2str(rall(2),3) ', p = ' num2str(p(2),3)])
% nexttile
% hold on
% plot(-1000:2000,mean(event_SDFs(matchedBySSDSeenTemp{1},:))*1000)
% plot(-1000:2000,mean(event_SDFs(matchedBySSDSeenTemp{2},:))*1000)
% xlim([-200 600])
% legend({'low','high'})
% nexttile
% hold on
% frs = sum(eventRasterstemp(:,1000:1200),2);
% 
% scatter(reg_tbl.RT(matchedBySSDSeenTemp{1}),frs(matchedBySSDSeenTemp{1}))
% scatter(reg_tbl.RT(matchedBySSDSeenTemp{2}),frs(matchedBySSDSeenTemp{2}))
[rerror,p] = corrcoef(reg_tbl.RT(matchedIndsTemp), frs(matchedIndsTemp));
% title(['r = ' num2str(rerror(2),3) ', p = ' num2str(p(2),3)])
reg_tbl.matched(matchedInds) =1;
errTable = reg_tbl_hilo(:,[1,2,6,7,9,11,12,24,25]);

glm_output.trial_errorHiLoSDF.rall = rall(2);
glm_output.trial_errorHiLoSDF.rerror = rerror(2);
glm_output.trial_errorHiLoSDF.CIs = errorCI;
glm_output.trial_errorHiLoSDF.diff = errorDiff;
glm_output.trial_errorHiLoSDF.table = actualDataHiLo;
glm_output.trial_errorHiLoSDF.ErrHi = ErrHi;
glm_output.trial_errorHiLoSDF.errLo = errLo;
glm_output.trial_errorHiLoSDF.corrHi = corrHi;
glm_output.trial_errorHiLoSDF.corrLo = corrLo;

%try again with RT matching

[matchedIndsTemp,matchedBySSDSeenTemp] = trlMatcherError(behData,goodTrlErrors,reg_tbl);
[matchedIndsCorr,matchedBySSDSeenCorr] = trlMatcherCorr(behData,goodTrlErrors,reg_tbl);
[~,event_SDFs,eventRasterstemp] = calc_FRs_matched_saccade(spk_data_in.Spikes,baseline_win,'saccade',...
    matchedBySSDSeenTemp,[-50, 350],fr_window,window_shift);
[~,event_SDFsCorr,eventRasterstempCorr] = calc_FRs_matched_saccade(spk_data_in.Spikes,baseline_win,'saccade',...
    matchedBySSDSeenCorr,[-50, 350],fr_window,window_shift);
eventRasterstempCorr(matchedIndsTemp,:) = eventRasterstemp(matchedIndsTemp,:);
event_SDFsCorr(matchedIndsTemp,:) = event_SDFs(matchedIndsTemp,:);
reg_tbl.matched_rasters3k(sort([matchedIndsTemp;matchedIndsCorr]),:) =...
    eventRasterstempCorr(sort([matchedIndsTemp;matchedIndsCorr]),:);
reg_tbl.sdfs3k(sort([matchedIndsTemp;matchedIndsCorr]),:) =...
    event_SDFsCorr(sort([matchedIndsTemp;matchedIndsCorr]),:);
reg_tbl_hiloRT = reg_tbl(sort([matchedIndsTemp;matchedIndsCorr]),:);

[errorCI, errorDiff, ~,actualDataHiLo] = ...
    calculatePermutationContinuity(reg_tbl_hiloRT.sdfs3k,...
    reg_tbl_hiloRT.trial_type==3 & reg_tbl_hiloRT.isH,....
    reg_tbl_hiloRT.trial_type==3 & ~reg_tbl_hiloRT.isH,20);


ErrHi = mean(reg_tbl_hiloRT.sdfs3k(reg_tbl_hiloRT.trial_type==3 & reg_tbl_hiloRT.isH,:));
errLo = mean(reg_tbl_hiloRT.sdfs3k(reg_tbl_hiloRT.trial_type==3 & ~reg_tbl_hiloRT.isH,:));
corrHi = mean(reg_tbl_hiloRT.sdfs3k(reg_tbl_hiloRT.trial_type==1 & reg_tbl_hiloRT.isH,:));
corrLo = mean(reg_tbl_hiloRT.sdfs3k(reg_tbl_hiloRT.trial_type==1 & ~reg_tbl_hiloRT.isH,:));

frs = sum(eventRasters(:,1000:1200),2);

[rall,p] = corrcoef(reg_tbl.RT(matchedInds), frs(matchedInds));

[rerror,p] = corrcoef(reg_tbl.RT(matchedIndsTemp), frs(matchedIndsTemp));
reg_tbl.matched(matchedInds) =1;
% errTable = reg_tbl_hilo(:,[1,2,6,7,9,11,12,24,25]);

glm_output.trial_errorHiLoSDFRT.rall = rall(2);
glm_output.trial_errorHiLoSDFRT.rerror = rerror(2);
glm_output.trial_errorHiLoSDFRT.CIs = errorCI;
glm_output.trial_errorHiLoSDFRT.diff = errorDiff;
glm_output.trial_errorHiLoSDFRT.table = actualDataHiLo;
glm_output.trial_errorHiLoSDFRT.ErrHi = ErrHi;
glm_output.trial_errorHiLoSDFRT.errLo = errLo;
glm_output.trial_errorHiLoSDFRT.corrHi = corrHi;
glm_output.trial_errorHiLoSDFRT.corrLo = corrLo;
%% corrected Error vs nostop
%take the stop signal aligned cancelled vs no stop, and subtract it from
%the stop SSRT aligned error, then realign on saccade and compare error
%to no stop

[~,event_SDFs,~] = calc_FRs_matched_saccade_sub_stop_signal(spk_data_in.Spikes,baseline_win,'saccade',...
    matchedBySSDSeen,[-50, 350],fr_window,window_shift,cancGoDiff,reg_tbl.RT,reg_tbl.SS,round(SSRT));

reg_tbl.sdfsCorr(matchedInds,:) = event_SDFs(matchedInds,:);
reg_tbl.matched_rasters(matchedInds,:) = eventRasters(matchedInds,1:2001);

reg_tbl_error = reg_tbl(matchedInds,:);

[errorCI, errorDiff, ~,actualDataErrorCorr] = ...
                calculatePermutationContinuity(reg_tbl_error.sdfsCorr,...
                reg_tbl_error.trial_type==1,reg_tbl_error.trial_type==3,contFiller);
noStopSDFErr = mean(reg_tbl_error.sdfsCorr(reg_tbl_error.trial_type==1,:));
nonCancSDFErr = mean(reg_tbl_error.sdfsCorr(reg_tbl_error.trial_type==3,:));

glm_output.trial_errorGoSDFCorr.CIs = errorCI;
glm_output.trial_errorGoSDFCorr.diff = errorDiff;
glm_output.trial_errorGoSDFCorr.table = actualDataErrorCorr;
glm_output.trial_errorGoSDFCorr.nostopSDF = noStopSDFErr;
glm_output.trial_errorGoSDFCorr.noCancSDF = nonCancSDFErr;


%calculate firing rates for the significant windows
%have to subtract the average activity for latency matched no stop
nostopFRs = [];

if ~isempty(actualData.modTimes)


    for ssd_i = find(cell2mat(cellfun(@(x) length(x)>4,matchedBySSD(:,1),'uni',0)))'
        for per = 1:length(actualData.mod1SDDur)
            tempNoStopFRs = sum(reg_tbl_trialtype.matched_rasters(...
                reg_tbl_trialtype.ssd == ssd_i & reg_tbl_trialtype.trial_type ==1,...
                actualData.modTimes{per}(1):actualData.modTimes{per}(2)),2)./ ...
                actualData.modDur{per}.* 1000;
            noStopFRsall{ssd_i,per} = tempNoStopFRs;
            nostopFRs(ssd_i,per) = mean(tempNoStopFRs);
            cancFR{ssd_i,per} = sum(reg_tbl_trialtype.matched_rasters(...
                reg_tbl_trialtype.ssd == ssd_i & reg_tbl_trialtype.trial_type ==2,...
                actualData.modTimes{per}(1):actualData.modTimes{per}(2)),2)./...
                actualData.modDur{per} .*1000;

        end
        tempFRs = ([cancFR{ssd_i,:}] - nostopFRs(ssd_i,:))./max(nostopFRs(ssd_i,:),.1);
        tempFRsNS = ([noStopFRsall{ssd_i,:}] - nostopFRs(ssd_i,:))./max(nostopFRs(ssd_i,:),.1);

        reg_tbl_trialtype.fr_diff(reg_tbl_trialtype.ssd == ssd_i & ...
            reg_tbl_trialtype.trial_type ==2,:) = tempFRs;
        reg_tbl_trialtype.fr_diff(reg_tbl_trialtype.ssd == ssd_i & ...
            reg_tbl_trialtype.trial_type ==1,:) = tempFRsNS;
        temp = mean([cancFR{ssd_i,:}] - nostopFRs(ssd_i,:))./max(nostopFRs(ssd_i,:),.1);
        glm_output.ssd.frDiffs{ssd_i,:} = temp;
    end
    nsFRs = reg_tbl_trialtype.fr_diff(reg_tbl_trialtype.trial_type ==1,:);
    cancFRs = reg_tbl_trialtype.fr_diff(reg_tbl_trialtype.trial_type ==2,:);

        glm_output.trial_typeSDF.AUC = calculateROCAUC(nsFRs,cancFRs);       

end

%calculate firing rates for the significant windows
%have to subtract the average activity for latency matched no stop
nostopFRs = [];
errorFR ={};
tempNoStopFRs= {};
if ~isempty(actualDataError.modTimes)


        for per = 1:length(actualDataError.mod1SDDur)
            tempNoStopFRs{per} = sum(reg_tbl_error.matched_rasters3k(...
                reg_tbl_error.trial_type ==1,...
                actualDataError.slopeTime{per}(1):actualDataError.slopeTime{per}(2)),2)./ ...
                actualDataError.slopeDur{per}.* 1000;
            nostopFRs(per) = mean(tempNoStopFRs{per});
            errorFR{1,per} = sum(reg_tbl_error.matched_rasters3k(...
                 reg_tbl_error.trial_type ==3,...
                actualDataError.slopeTime{per}(1):actualDataError.slopeTime{per}(2)),2)./...
                actualDataError.slopeDur{per} .*1000;
            %2SD
 tempNoStopFRs2sd{per} = sum(reg_tbl_error.matched_rasters3k(...
                reg_tbl_error.trial_type ==1,...
                actualDataError.modTimes{per}(1):actualDataError.modTimes{per}(2)),2)./ ...
                actualDataError.modDur{per}.* 1000;
            errorFR2sd{1,per} = sum(reg_tbl_error.matched_rasters3k(...
                 reg_tbl_error.trial_type ==3,...
                actualDataError.modTimes{per}(1):actualDataError.modTimes{per}(2)),2)./...
                actualDataError.modDur{per} .*1000;
            hiErr{per} = sum(reg_tbl_hilo.matched_rasters3k(...
                reg_tbl_hilo.trial_type ==3 & reg_tbl_hilo.isH==1,...
                actualDataError.modTimes{per}(1):actualDataError.modTimes{per}(2)),2)./ ...
                actualDataError.slopeDur{per}.* 1000;
            loErr{per} = sum(reg_tbl_hilo.matched_rasters3k(...
                reg_tbl_hilo.trial_type ==3 & reg_tbl_hilo.isL==1,...
                actualDataError.modTimes{per}(1):actualDataError.modTimes{per}(2)),2)./ ...
                actualDataError.slopeDur{per}.* 1000;
            hiCorr{per} = sum(reg_tbl_hilo.matched_rasters3k(...
                reg_tbl_hilo.trial_type ==1 & reg_tbl_hilo.isH==1,...
                actualDataError.modTimes{per}(1):actualDataError.modTimes{per}(2)),2)./ ...
                actualDataError.slopeDur{per}.* 1000;
            loCorr{per} = sum(reg_tbl_hilo.matched_rasters3k(...
                reg_tbl_hilo.trial_type ==1 & reg_tbl_hilo.isL==1,...
                actualDataError.modTimes{per}(1):actualDataError.modTimes{per}(2)),2)./ ...
                actualDataError.slopeDur{per}.* 1000;

            %narrow
             tempNoStopFRsNarr{per} = sum(reg_tbl_error.matched_rasters3k(...
                reg_tbl_error.trial_type ==1,...
                actualDataError.narrTime{per}(1):actualDataError.narrTime{per}(2)),2)./ ...
                actualDataError.narrDur{per}.* 1000;
            errorFRNarr{1,per} = sum(reg_tbl_error.matched_rasters3k(...
                 reg_tbl_error.trial_type ==3,...
                actualDataError.narrTime{per}(1):actualDataError.narrTime{per}(2)),2)./...
                actualDataError.narrDur{per} .*1000;
        end

        glm_output.trial_errorGoSDF.AUC = calculateROCAUC(cell2mat(tempNoStopFRs),cell2mat(errorFR));       
        glm_output.trial_errorGoSDF.AUC2sd = calculateROCAUC(cell2mat(tempNoStopFRs2sd),cell2mat(errorFR2sd));
        glm_output.trial_errorGoSDF.AUCNarr = calculateROCAUC(cell2mat(tempNoStopFRsNarr),cell2mat(errorFRNarr));
        reg_tbl_error.fr_diff( reg_tbl_error.trial_type ==3,:) = [errorFR{1,:}] - nostopFRs(1,:);
        glm_output.nonCanc.frDiffs{1,:} = mean([errorFR{1,:}] - nostopFRs(1,:));
        reg_tbl_error.ssdVal(reg_tbl_error.ssd==ssd_i,:) = inh_pnc(2,ssd_i);
        reg_tbl_error.pnc(reg_tbl_error.ssd==ssd_i,:) = inh_pnc(1,ssd_i);
        HiLotable = table();

        
        HiLotable.errHiLop = cell2mat(arrayfun(@(x) ranksum(hiErr{x},loErr{x}),1:per,'uni',0)');
        HiLotable.CorrHiLop = cell2mat(arrayfun(@(x) ranksum(hiCorr{x},loCorr{x}),1:per,'uni',0)');
        HiLotable.errHiLoD = cell2mat(arrayfun(@(x) cohens_d(hiErr{x},loErr{x}),1:per,'uni',0)');
        HiLotable.CorrHiLoD = cell2mat(arrayfun(@(x) cohens_d(hiCorr{x},loCorr{x}),1:per,'uni',0)');
        HiLotable.ErrValInd = cell2mat(arrayfun(@(x) (mean(hiErr{x}) -mean(loErr{x}))/...
            (mean(hiErr{x}) + mean(loErr{x})),1:per,'uni',0)');
        HiLotable.CorrValInd = cell2mat(arrayfun(@(x) (mean(hiCorr{x}) -mean(loCorr{x}))/...
            (mean(hiCorr{x}) + mean(loCorr{x})),1:per,'uni',0)');
        HiLotable.errN = cellfun(@length, hiErr)';
        HiLotable.corrN = cellfun(@length, hiCorr)';
        glm_output.trial_errorGoSDF.hiloTable = HiLotable;
else
     glm_output.trial_errorGoSDF.AUC =[];
     glm_output.trial_errorGoSDF.AUC2sd =[];
     glm_output.trial_errorGoSDF.AUCNarr =[];
end
%% Run GLM: stop_signal delay

    sigPermTimes = length(actualData.modDur);

glm_output.ssd.sig_times = zeros(1,sigPermTimes); 
glm_output.ssd.beta_weights = zeros(1,sigPermTimes); 
glm_output.ssd.permPos = zeros(1,sigPermTimes);
glm_output.ssd.permSig = zeros(1,sigPermTimes);
glm_output.ssd.pvals = zeros(1,sigPermTimes);
u_t_mdl = [];

reg_tbl_stopping = reg_tbl_trialtype(reg_tbl_trialtype.is_canc==1,:);
%check that there are enough trials of eac SSD
usedSSDs = unique(reg_tbl_stopping.ssd);
for ssd_i = 1:length(usedSSDs)
    if height(reg_tbl_stopping(reg_tbl_stopping.ssd==usedSSDs(ssd_i),:))<5
        %remove the SSDs that have too few trials
        reg_tbl_stopping(reg_tbl_stopping.ssd==usedSSDs(ssd_i),:) = [];
        glm_output.ssd.cancSDFs{ssd_i} = [];
        glm_output.ssd.nsSDFs{ssd_i} = [];
    else
        glm_output.ssd.cancSDFs{ssd_i} = mean(reg_tbl_stopping.sdfs(...
            reg_tbl_stopping.ssd==usedSSDs(ssd_i) & reg_tbl_stopping.is_canc ==1,:));
        glm_output.ssd.nsSDFs{ssd_i} = mean(reg_tbl_trialtype.sdfs(...
            reg_tbl_trialtype.ssd==usedSSDs(ssd_i) & reg_tbl_trialtype.is_canc ==0,:));
    end
end

    sigPermTimes = length(actualData.modDur);


% For each averaged time point
for timepoint_i = 1:sigPermTimes
    
    % Input the timepoint specific firing times
    reg_tbl_stopping.firing_rate = reg_tbl_stopping.fr_diff(:,timepoint_i); 
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ ssdVal + trl_numZ');
    
    % GLM output ---------------------------------   
    % - Stop-signal delay
    glm_output.ssd.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % ssd
    glm_output.ssd.pvals(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2); % ssd
    glm_output.ssd.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.Estimate(2); % ssd
    
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
glm_output.pnc.sig_times = zeros(1,sigPermTimes); 
glm_output.pnc.beta_weights = zeros(1,sigPermTimes); 
glm_output.pnc.permPos = zeros(1,sigPermTimes);
glm_output.pnc.permSig = zeros(1,sigPermTimes);
glm_output.pnc.pvals = zeros(1,sigPermTimes);
u_t_mdl = [];


% For each averaged time point
for timepoint_i = 1:sigPermTimes
    
    % Input the timepoint specific firing times
    reg_tbl_stopping.firing_rate = reg_tbl_stopping.fr_diff(:,timepoint_i); 
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ pnc + trl_numZ');
    
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

  sigPermTimes = length(actualData.modDur);

glm_output.surp.sig_times = zeros(1,sigPermTimes); 
glm_output.surp.beta_weights = zeros(1,sigPermTimes); 
glm_output.surp.permPos = zeros(1,sigPermTimes);
glm_output.surp.permSig = zeros(1,sigPermTimes);
glm_output.surp.pvals = zeros(1,sigPermTimes);
u_t_mdl = [];

reg_tbl_stopping = reg_tbl_trialtype(reg_tbl_trialtype.is_canc==1,:);
%check that there are enough trials of eac SSD
usedSSDs = unique(reg_tbl_stopping.ssd);
for ssd_i = 1:length(usedSSDs)
    if height(reg_tbl_stopping(reg_tbl_stopping.ssd==usedSSDs(ssd_i),:))<5
        reg_tbl_stopping(reg_tbl_stopping.ssd==usedSSDs(ssd_i),:) = [];
    end
end

    sigPermTimes = length(actualData.modDur);


% For each averaged time point
for timepoint_i = 1:sigPermTimes
    
    % Input the timepoint specific firing times
    reg_tbl_stopping.firing_rate = reg_tbl_stopping.fr_diff(:,timepoint_i); 
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ surprise + trl_numZ');
    
    % GLM output ---------------------------------   
    % - Stop-signal delay
    glm_output.surp.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % ssd
    glm_output.surp.pvals(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2); % ssd
    glm_output.surp.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.Estimate(2); % ssd
    
%     % - Value
%     glm_output.ssd.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % value
%     glm_output.ssd.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % value
%         
%     % - SSD x Value 
%     glm_output.ssd.sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(5) < .01; % ssd x value
%     glm_output.ssd.beta_weights(3,timepoint_i) = u_t_mdl.Coefficients.tStat(5); % % ssd x value
    
    % - Trial number
    glm_output.surp.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.surp.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.Estimate(3); % trial_number

    
end
%% determine if there is a difference between cancelled and noncancelled
glm_output.nonCanc.sig_times = []; 
glm_output.nonCanc.beta_weights = []; 
u_t_mdl = [];

%% seen vs unseen
glm_output.seen.sig_times = []; 
glm_output.seen.beta_weights = []; 
u_t_mdl = []; 
trlInds = reg_tbl.trial_number;
%% Run GLM: stop_signal delay

    sigPermTimes = length(actualDataError.modDur);

glm_output.ssderror.sig_times = zeros(1,sigPermTimes); 
glm_output.ssderror.beta_weights = zeros(1,sigPermTimes); 
glm_output.ssderror.permPos = zeros(1,sigPermTimes);
glm_output.ssderror.permSig = zeros(1,sigPermTimes);
glm_output.ssderror.pvals = zeros(1,sigPermTimes);
u_t_mdl = [];

reg_tbl_error = reg_tbl_error(reg_tbl_error.trial_type ==3,:);
%check that there are enough trials of eac SSD
usedSSDs = unique(reg_tbl_error.ssd);
for ssd_i = 1:length(usedSSDs)
    if height(reg_tbl_error(reg_tbl_error.ssd==usedSSDs(ssd_i),:))<5
        reg_tbl_error(reg_tbl_error.ssd==usedSSDs(ssd_i),:) = [];
    end
end



% For each averaged time point
for timepoint_i = 1:sigPermTimes
    
    % Input the timepoint specific firing times
    reg_tbl_error.firing_rate = reg_tbl_error.fr_diff(:,timepoint_i); 
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_error,'firing_rate ~ ssdVal + trl_numZ');
    
    % GLM output ---------------------------------   
    % - Stop-signal delay
    glm_output.ssderror.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < regAlpha; % ssd
    glm_output.ssderror.pvals(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2); % ssd
    glm_output.ssderror.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.Estimate(2); % ssd
    
%     % - Value
%     glm_output.ssd.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % value
%     glm_output.ssd.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % value
%         
%     % - SSD x Value
%     glm_output.ssd.sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(5) < .01; % ssd x value
%     glm_output.ssd.beta_weights(3,timepoint_i) = u_t_mdl.Coefficients.tStat(5); % % ssd x value
    
    % - Trial number
    glm_output.ssderror.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.ssderror.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.Estimate(3); % trial_number

    
end

%% Determine periods of significance
%if seen isn't significant, we don't care about SSD being significant
% glm_output.ssd.sig_times(1,glm_output.seen.sig_times(1,:) ==0)=0;
% 
%if conflict isn't significant, we don't care about SSD being significant
% glm_output.ssd.sig_times(1,glm_output.trial_type.sig_times(1,:) ==0)=0;
% glm_output.ssd.sig_timePos = glm_output.ssd.sig_times & (glm_output.ssd.beta_weights>0);
% glm_output.pnc.sig_times(1,glm_output.trial_type.sig_times(1,:) ==0)=0;
% glm_output.pnc.sig_timePos = glm_output.pnc.sig_times & (glm_output.pnc.beta_weights>0);

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

% 
% signal_detect_length = 50;
% signal_detect_wins = signal_detect_length/window_shift;
% 
% % now ask whether this unit was significant with significance defined as at least
% % 100ms of selectivity following SSD
% 
% % Trial-type
% [canc_sig_start, canc_sig_len, ~] = ZeroOnesCount(glm_output.trial_type.sig_times(1,:)); % choice direction
% 
% % Stop-signal delay
% [ssd_sig_start, ssd_sig_len, ~] = ZeroOnesCount(glm_output.ssd.sig_timePos(1,:) ); % choice direction
% 
% % inhibition function
% [pnc_sig_start, pnc_sig_len, ~] = ZeroOnesCount(glm_output.pnc.sig_timePos(1,:) ); % choice direction
% 
% % nonCanc
% [~, noncanc_sig_len, ~] = ZeroOnesCount(glm_output.nonCanc.sig_times(1,:)); % choice direction
% 
% % Seen
% [~, seen_sig_len, ~] = ZeroOnesCount(glm_output.seen.sig_times(1,:)); % choice direction
% 
% % SSD perm
% [~, ssdP_sig_len, ~] = ZeroOnesCount(glm_output.ssd.permSig(1,:)); % choice direction
% 
% 
% encoding_flag = []; encoding_beta = [];
% 

encoding_flag(1,1) = ~isempty(glm_output.trial_typeSDF.table.modDur);
encoding_flag(1,2) = any(glm_output.ssd.sig_times(1,:));
encoding_flag(1,3) = any(glm_output.pnc.sig_times(1,:));
encoding_flag(1,4) = any(glm_output.surp.sig_times(1,:));
encoding_flag(1,5) = ~isempty(glm_output.trial_canNCSDF.table.modDur);
encoding_flag(1,6) = ~isempty(glm_output.trial_errorGoSDF.table.modDur);
encoding_flag(1,7) = any(glm_output.ssderror.sig_times(1,:));
% encoding_flag(1,5) = any(seen_sig_len >= signal_detect_wins);
% encoding_flag(1,6) = any(ssdP_sig_len >= signal_detect_wins);
% sig_bins = cell(1,3);
% if encoding_flag(1,1)
%     cancInd = find(canc_sig_len>=5,1);
%     sig_bins{1,1} = canc_sig_start(cancInd):canc_sig_start(cancInd)+canc_sig_len(cancInd)-1;
%     FRs = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
%         [-60+(10*sig_bins{1}(1)) -60 + 100 + (10*sig_bins{1}(end))],canc_sig_len(cancInd)*10+90,window_shift);
%     frs.canc{1} = mean(FRs(reg_tbl_trialtype.trial_number(reg_tbl_trialtype.trial_type==1)));
%     frs.canc{2}= mean(FRs(reg_tbl_trialtype.trial_number(reg_tbl_trialtype.trial_type==2)));
% end
% pvals = struct();
% if encoding_flag(1,2) 
%     ssdInd = find(ssd_sig_len>=5,1);
%     sig_bins{1,2} = ssd_sig_start(ssdInd):ssd_sig_start(ssdInd)+ssd_sig_len(ssdInd)-1;
%     FRs = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
%         [-60+(10*sig_bins{2}(1)) -60 + 100 + (10*sig_bins{2}(end))],ssd_sig_len(ssdInd)*10+90,window_shift);
%     ssds = unique(reg_tbl_stopping.ssd);
%     for ii = 1:length(ssds)
%         latMatched = mean(FRs(reg_tbl_trialtype.trial_number(...
%             reg_tbl_trialtype.ssd == ssds(ii) & reg_tbl_trialtype.trial_type ==1)));
%         frs.ssd{ii} = mean(FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))));
%         tempfrs{ii} = FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))) - latMatched;
%         tempssd{ii}  = ones(length(tempfrs{ii}),1)*inh_pnc(2,ssds(ii));
%         tempPNC{ii} = ones(length(tempfrs{ii}),1)*inh_pnc(1,ssds(ii));
%         frs.ssdDiff{ii} = frs.ssd{ii}-latMatched;
%     end
%     ssdfitssd = fitlm(vertcat(tempssd{:}), vertcat(tempfrs{:}));
%     ssdfitssdp = ssdfitssd.anova.pValue;
%     ssdfitpnc = fitlm(vertcat(tempPNC{:}), vertcat(tempfrs{:}));
%     ssdfitpncp = ssdfitpnc.anova.pValue;
%     pvals.SSD = table(ssdfitssdp,ssdfitpncp);
% end
% if encoding_flag(1,3)
%     pncInd = find(pnc_sig_len>=5,1);
%     sig_bins{1,3} = pnc_sig_start(pncInd):pnc_sig_start(pncInd)+pnc_sig_len(pncInd)-1;
% 
%     FRs = calc_FRs_and_zscore(spk_data_in.Spikes,baseline_win,event_alignment,...
%         [-60+(10*sig_bins{3}(1)) -60 + 100 + (10*sig_bins{3}(end))],pnc_sig_len(pncInd)*10+90,window_shift);
%     ssds = unique(reg_tbl_stopping.ssd);
%     for ii = 1:length(ssds)
%          latMatched = mean(FRs(reg_tbl_trialtype.trial_number(...
%             reg_tbl_trialtype.ssd == ssds(ii) & reg_tbl_trialtype.trial_type ==1)));
% 
%         frs.pnc{ii} = mean(FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))));
%         frs.pncDiff{ii} = frs.pnc{ii} - latMatched;
%          tempfrs{ii} = FRs(reg_tbl_stopping.trial_number(reg_tbl_stopping.ssd==ssds(ii))) - latMatched;
%         tempssd{ii}  = ones(length(tempfrs{ii}),1)*inh_pnc(2,ssds(ii));
%         tempPNC{ii} = ones(length(tempfrs{ii}),1)*inh_pnc(1,ssds(ii));
%     end
%     pncfitssd = fitlm(vertcat(tempssd{:}), vertcat(tempfrs{:}));
%     pncfitssdp = pncfitssd.anova.pValue;
%     pncfitpnc = fitlm(vertcat(tempPNC{:}), vertcat(tempfrs{:}));
%     pncfitpncp = pncfitpnc.anova.pValue;
%     pvals.PNC = table(pncfitssdp,pncfitpncp);
% end
% % Get average beta-weights
% encoding_beta(1,1) = nanmean(glm_output.trial_type.beta_weights(1,sig_bins{1}));
% encoding_beta(1,2) = nanmean(glm_output.ssd.beta_weights(1,sig_bins{2}));
% encoding_beta(1,3) = nanmean(glm_output.pnc.beta_weights(1,sig_bins{3}));
% encoding_beta(1,4) = nanmean(glm_output.nonCanc.beta_weights(1,...
%     glm_output.nonCanc.sig_times(1,:)==1));
% encoding_beta(1,5) = nanmean(glm_output.seen.beta_weights(1,...
%     glm_output.seen.sig_times(1,:)==1));



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
