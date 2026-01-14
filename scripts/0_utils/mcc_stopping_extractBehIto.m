function [behavior,mcc_map_info] = mcc_stopping_extractBehIto(dirs,dataFiles_beh)
a = 0;
counter = 1;
mcc_map_info = table();

% Looping through each of the individual data files
for dataFileIdx = 1:length(dataFiles_beh)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_beh{dataFileIdx},dataFileIdx,length(dataFiles_beh))
    
    % We then get the (behavior) filename of the record of interest
    behFilename = [dataFiles_beh{dataFileIdx} '-beh'];
    % and load it into the workspace
    import_data = struct(); import_data = load([dirs dataFiles_beh{dataFileIdx}]);
    % trialsTbl = dataset2table(import_data.trialData)
    trialsTbl = import_data.trialData;
        sessData = import_data.sessionData;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Once we have the behavior in matlab, we can then look to extract
    % relevant task/session behavior.
    
    sessionName = {behFilename}; % Get the session name
    [ttx, ttx_history, trialEventTimes] =... % Index of trials and event timings
        beh_getTrialsIto(trialsTbl);
    [stopSignalBeh, ~] = beh_getStoppingInfoIto... % Stopping behavior
        (trialsTbl,trialEventTimes,ttx);
    if min(stopSignalBeh.inh_weibull.y)<=0.2 && max(stopSignalBeh.inh_weibull.y)>=0.8...
            && stopSignalBeh.ssrt.integrationWeighted > 60

        trialEventTimes.ssrt = trialEventTimes.stopSignal_artifical + stopSignalBeh.ssrt.integrationWeighted;
        [ttm] = beh_getMatchedTrials(stopSignalBeh,ttx, trialEventTimes); % Trial matching indices
        [ttseen] = beh_getSeenTrials(stopSignalBeh,ttx, trialEventTimes); % Trial matching indices

        % After extracting the individual behavioral variable, we then collapse
        % it into one structure for the given session.
        behavior(counter) = struct('sessionName',sessionName,'ttx',ttx,'trialEventTimes',trialEventTimes,...
            'stopSignalBeh',stopSignalBeh,'ttm',ttm,'ttx_history',ttx_history,'ttseen',ttseen);
        temp_mcc_map = processNoHoNeuralData(trialsTbl,trialEventTimes,sessData);
    

    mcc_map_info = [mcc_map_info;temp_mcc_map];
        counter = counter+1;
    end
end

behavior = struct2table(behavior);


end

