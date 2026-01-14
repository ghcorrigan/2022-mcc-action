function [mcc_map_info,bankInfo] = mcc_dv_mapChNewBank(dajo_datamap_curated, dirs)
% Get session information and indices
n_sessions = size(dajo_datamap_curated,1);
mcc_map_info = [];
bankInfo = struct();

for session_i = 1:n_sessions
    neuralFilename = dajo_datamap_curated.dataFilename{session_i};
    % behFilename = data_findBehFile(neuralFilename);
    % beh_index = util_find_beh_index(behavior,behFilename);
    logInfo(session_i,:) = util_getLogInfo(neuralFilename);
    
    fprintf(['Extracting data for ' neuralFilename ': session %i of %i     \n'],...
        session_i, n_sessions);
    
    acs_ch_mapping = util_getACCchannels(logInfo,session_i);
    
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
    probeName = probeName{1}(1:end-4);
    bankInfo.(probeName) =table();
    
    %check whether there are any units on the bank channel
    if any(map_info.depth==0) && any(map_info.depth)
        spacing = abs(diff(acs_ch_mapping.um_depth_acs(1:2)));
       %check which direction has the larger gap, if there is one
       depths = map_info.depth;
       zeroChans = map_info.depth==0;
       %if the first unit is in the bank, assign it to the ventral bank,
       %and if it is the last unit, assign it to the dorsal bank, but
       %otherwise measure the distance to the next unit in each direction
       %and if there is a closer bank, assign it to that bank
       if zeroChans(1)==1
           preGap=1000;
       else
       preGap = abs(depths(find(zeroChans,1)-1));
       end
       if zeroChans(end)==1
           postGap = 1000;
       else
       postGap = abs(depths(find(zeroChans,1,'last')+1));
       end
       if preGap==postGap
           %do nothing
       elseif preGap>postGap
           % change the bank channel to the channel above the current one
           map_info.depth = map_info.depth + spacing;

       elseif postGap>preGap
           % change the bank channel to the channel below the current one
           map_info.depth = map_info.depth - spacing;

       end

    end
    map_info.area(map_info.depth<0) = {'dMCC'};
    map_info.area(map_info.depth>0) = {'vMCC'};
    bankInfo.(probeName).Units = map_info.unit;
    bankInfo.(probeName).PUnits = map_info.probeUnit;
    bankInfo.(probeName).dACC = map_info.depth<0;
    bankInfo.(probeName).vACC = map_info.depth>0;
    bankInfo.(probeName).ambiguous = map_info.depth<151 & map_info.depth>-151;


    mcc_map_info = [mcc_map_info; map_info];
end

% Clear up unknown areas and remove noise clusters
mcc_map_info(strcmp(mcc_map_info.area,'?'),:) = [];
mcc_map_info(mcc_map_info.noise == 1,:) = [];

% dataFiles_beh = unique(behavior.sessionName);
% 
% % Report numbers
% fprintf('Nsessions, Monkey Da: %i | Nsessions, Monkey Jo: %i    \n',...
%     sum(contains(dataFiles_beh,'dar')), sum(contains(dataFiles_beh,'jou')))
fprintf('Nneurons [all], Monkey Da: %i | Nneurons [all], Monkey Jo: %i    \n',...
    sum(contains(mcc_map_info.monkey,'dar')), sum(contains(mcc_map_info.monkey,'jou')))
fprintf('Nneurons [dMCC], Monkey Da: %i | Nneurons [dMCC], Monkey Jo: %i    \n',...
    sum(contains(mcc_map_info.monkey,'dar') & contains(mcc_map_info.area,'dMCC')),...
    sum(contains(mcc_map_info.monkey,'jou') & contains(mcc_map_info.area,'dMCC')))
fprintf('Nneurons [vMCC], Monkey Da: %i | Nneurons [vMCC], Monkey Jo: %i    \n',...
    sum(contains(mcc_map_info.monkey,'dar') & contains(mcc_map_info.area,'vMCC')),...
    sum(contains(mcc_map_info.monkey,'jou') & contains(mcc_map_info.area,'vMCC')))
fprintf('Nneurons [ACS], Monkey Da: %i | Nneurons [ACS], Monkey Jo: %i    \n',...
    sum(contains(mcc_map_info.monkey,'dar') & contains(mcc_map_info.area,'ACS')),...
    sum(contains(mcc_map_info.monkey,'jou') & contains(mcc_map_info.area,'ACS')))
