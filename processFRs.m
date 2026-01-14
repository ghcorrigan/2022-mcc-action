for neuron = 1:length(cluster_neurons)
    neuron_i = cluster_neurons(neuron);
   neuralFilename = mcc_map_info.session{neuron_i};%... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
    warning off
    [frs.(['neuron_' num2str(neuron_i)])] = comp_FR_rates(neuron_i, mcc_map_info, dataFiles_beh, dirs,behavior(behaviorIdx,:));
    
    
end