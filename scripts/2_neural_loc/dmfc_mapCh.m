function dmfc_map_info = dmfc_mapCh(dajo_datamap_curated, behavior, dirs)
% Get session information and indices
n_sessions = size(dajo_datamap_curated,1);
dmfc_map_info = [];

for session_i = 1:n_sessions
    neuralFilename = dajo_datamap_curated.dataFilename{session_i};
    behFilename = data_findBehFile(neuralFilename);
    beh_index = util_find_beh_index(behavior,behFilename);
    logInfo(session_i,:) = util_getLogInfo(neuralFilename);
    
    fprintf(['Extracting data for ' neuralFilename ': session %i of %i     \n'],...
        session_i, n_sessions);
    
    acs_ch_mapping = util_getDMFCchannels(logInfo,session_i);
    
    spk_data = load(fullfile(dirs.data,[neuralFilename '-spk.mat']));
    
    n_neurons= size(dajo_datamap_curated.spkInfo(session_i,1).unitInfo,1);
    map_info = table();
    
    for neuron_i = 1:n_neurons
        site = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.site(neuron_i);
        session = dajo_datamap_curated.dataFilename(session_i);
        monkey = dajo_datamap_curated.monkey(session_i);
        unit = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.unitDSP(neuron_i);
        depth = acs_ch_mapping.um_depth_acs(site);
        area = acs_ch_mapping.area(site);
        ap = logInfo.ap_stereo(session_i);
        ml = logInfo.ml_stereo(session_i);
        mua = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.flag_mua(neuron_i);
        noise = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.flag_noise(neuron_i);

        spk_width = util_getSpkWidth(spk_data,dajo_datamap_curated.spkInfo(session_i,1).unitInfo.unitWAV{neuron_i});
        probeName = dajo_datamap_curated.probes(session_i);
        probeNum = str2double(probeName{1}(7:end-4));
        probeUnit = neuron_i;

        map_info(neuron_i,:) = table(session,monkey,unit,site,depth,area,ap,ml,mua,noise,spk_width,probeName,probeNum,probeUnit);
    end
    
    dmfc_map_info = [dmfc_map_info; map_info];
end

% Clear up unknown areas and remove noise clusters
% dmfc_map_info(strcmp(dmfc_map_info.area,'?'),:) = [];
dmfc_map_info(dmfc_map_info.noise == 1,:) = [];

dmfc_map_info.area = repmat({'DMFC'},size(dmfc_map_info,1),1);

dataFiles_beh = unique(behavior.sessionName);

% Report numbers
fprintf('Nsessions, Monkey Da: %i | Nsessions, Monkey Jo: %i    \n',...
    sum(contains(dataFiles_beh,'dar')), sum(contains(dataFiles_beh,'jou')))
fprintf('Nneurons [all], Monkey Da: %i | Nneurons [all], Monkey Jo: %i    \n',...
    sum(contains(dmfc_map_info.monkey,'dar')), sum(contains(dmfc_map_info.monkey,'jou')))
