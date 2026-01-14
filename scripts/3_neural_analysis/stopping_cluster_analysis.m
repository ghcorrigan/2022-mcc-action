
%% Setup & parameterize clustering
input_neurons = [];
input_neurons = find(mcc_analysis_table.glm_trial == 1);

canc_sdf_in = []; nostop_sdf_in = []; noncanc_sdf_in = [];

for neuron_i = 1:length(input_neurons)
    neuron_j = input_neurons(neuron_i);
    canc_sdf_in(neuron_i,:) = nanmean(mcc_analysis_table.sdf_canceled{neuron_j});
    nostop_sdf_in(neuron_i,:) = nanmean(mcc_analysis_table.sdf_nostop{neuron_j});
    noncanc_sdf_in(neuron_i,:) = nanmean(mcc_analysis_table.sdf_noncanc{neuron_j});

    diff_sdf_in(neuron_i,:) = canc_sdf_in(neuron_i,:)-nostop_sdf_in(neuron_i,:);
    diff_sdf_in(neuron_i,:) = diff_sdf_in(neuron_i,:)./max(diff_sdf_in(neuron_i,1000+[0:600]));
end

inputSDF = [];
inputSDF = {diff_sdf_in};

sdfTimes = {[-1000:1000]};
sdfEpoch = {[0:600]};

colorMapping = [1,1];

%% Clustering algorithm
[sortIDs,idxDist, raw, respSumStruct, rawLink,idK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0.5);

%% Cluster identification
% Refine number of clusters (iterative with below dendrogram)
myK = 10; % idK = based on consensusCluster raw output

%% Identify/assign neurons to clusters
nClusters_manual = myK; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end

%% Plot dendrogram
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,myK))
set(gca,'YTick',[],'YLim',[1 length(rawLink)]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
imagesc(raw(outPerm,outPerm));
colormap(flipud(summer));
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])

%% Plot cluster populations
% Generate a quick population sdf of each cluster
for cluster_i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    example_neuron_subfig = subplot(1,1,1); hold on
    plot(sdfTimes{1},nanmean(canc_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{1},nanmean(nostop_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
%     plot(sdfTimes{1},nanmean(noncanc_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.noncanc]);
    vline([0], 'k--'); xlim([-250 750])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
end


%% Merge clusters
% Merge: 4, 8 (transient)

% Keep: 3(suppressed), 6 (facilitated, sustained), 7 (transient
% suppressed), 9 (oscil).

% Remove: 1, 2, 5, 10

cluster_merge_idx = {[4,8],3,6,7,9};

for cluster_merge_i = 1:length(cluster_merge_idx)
    cluster_neuron_id{cluster_merge_i} = clusterNeurons{cluster_merge_idx{cluster_merge_i}};
    cluster_neuron_id{cluster_merge_i} = sort(cluster_neuron_id{cluster_merge_i});
end

%% Check merge population sdf
% Generate a quick population sdf of each cluster
for cluster_i = 1:length(cluster_neuron_id)
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    example_neuron_subfig = subplot(1,1,1); hold on
    plot(sdfTimes{1},nanmean(canc_sdf_in(cluster_neuron_id{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{1},nanmean(nostop_sdf_in(cluster_neuron_id{cluster_i},:),1), 'color', [colors.nostop]);
    plot(sdfTimes{1},nanmean(noncanc_sdf_in(cluster_neuron_id{cluster_i},:),1), 'color', [colors.noncanc]);
    vline([0], 'k--'); xlim([-250 750])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(cluster_neuron_id{cluster_i}))])
end

%% Plot individual cluster neurons
% Print each neuron within a cluster
for cluster_i = 1:length(cluster_neuron_id)
    
    input_neurons_cluster = [];
    input_neurons_cluster = cluster_neuron_id{cluster_i};
    
    n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
    n_batches = ceil(size(input_neurons_cluster,1)/n_plot_sheet);
    
    neuron_i = 0;
    for page_i = 1:n_batches
        fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
        
        for plot_i = 1:n_plot_sheet
            neuron_i = neuron_i+1;
            try
                neuron_j = input_neurons_cluster(neuron_i);
                
                neuralFilename = mcc_map_info.session{neuron_j};
                neuronLabel = mcc_map_info.unit{neuron_j};
                
                %... and find the corresponding behavior file index
                behFilename = data_findBehFile(neuralFilename);
                behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
                
                
                subplot(n_plot_x, n_plot_y, plot_i); hold on
                plot(-1000:2000, canc_sdf_in(neuron_j,:),'color',colors.canceled)
                plot(-1000:2000, nostop_sdf_in(neuron_j,:),'color',colors.nostop)
                xlim([-250 750]); vline(0,'k-')
                xlabel('Time from stop-signal (ms)')
                ylabel('Firing rate (spks/sec)')
                title(['Cluster: ' int2str(cluster_i) ' | Neuron: ' int2str(neuron_j) ])
                vline(behavior.stopSignalBeh(behaviorIdx).ssrt.integrationWeighted,'k--')
            catch
                continue
            end
            
        end
        
        filename = fullfile(dirs.root,'results','cluster',['cluster_canc_' int2str(cluster_i) '_pg' int2str(page_i) '.pdf']);
        set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
        print(fig_out,filename,'-dpdf') % then print it
        close(fig_out)
    end
end

%% Heatmap: Cluster activity
% Generate a quick population sdf of each cluster
figure('Renderer', 'painters', 'Position', [100 100 750 400]);hold on

n_clusters = length(cluster_neuron_id);

for cluster_i = 1:n_clusters
    sdf_in_diff = []; sdf_in_diff = diff_sdf_in(cluster_neuron_id{cluster_i},:);
    
    temp_a = []; peak_latency = []; neuron_order = [];
    temp_a = sdf_in_diff(:,1000+[0:600]);
    
    for neuron_i = 1:size(temp_a,1)
        peak_latency(neuron_i,1) = find(temp_a(neuron_i,:) == max(temp_a(neuron_i,:)),1);
    end
    [~,neuron_order] = sort(peak_latency,'ascend'); % sorted from max to min
    
    
    subplot(3,n_clusters,[cluster_i]); hold on
    plot(sdfTimes{1},nanmean(canc_sdf_in(cluster_neuron_id{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{1},nanmean(nostop_sdf_in(cluster_neuron_id{cluster_i},:),1), 'color', [colors.nostop]);
    vline([0], 'k--'); xlim([-250 750])
    
    subplot(3,n_clusters,[cluster_i+n_clusters cluster_i+n_clusters*2]); hold on
    imagesc('XData',[-1000:2000],'YData',[1:size(sdf_in_diff,1)],'CData',sdf_in_diff(neuron_order,:))
    colormap('viridis')
    vline([0], 'k--'); xlim([-200 600]); ylim([1-0.5 max(cellfun(@length,cluster_neuron_id))+0.5]); caxis([0 1])
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(cluster_neuron_id{cluster_i}))])
    set(gca,'YDir','Reverse')
    
end

%% Figure: Raster, example SDF, population SDF

% Example neuron idx ----------------------------
% Cluster 1: 1
% Cluster 2: 51
% Cluster 3: 73/166/246
% Cluster 4: 30
% Cluster 5: 32/105/199

example_neuron_idx = [1, 51, 166, 30, 105];
example_neuron_ylim = {[0 15],[0 15],[5 25],[0 10],[5 15]};
for cluster_i = 1:5
    input_neuron_idx = example_neuron_idx(cluster_i);
    raw_neuron_idx = input_neurons(input_neuron_idx);
    
    main_sdf_fig_out = [];
    main_sdf_fig_out = main_sdf_figure(dirs,colors,...
        mcc_analysis_table,raw_neuron_idx,behavior,dataFiles_beh,...
        clusterNeurons, cluster_i, canc_sdf_in, nostop_sdf_in);
    
    main_sdf_fig_out(2,1).axe_property('YLim',example_neuron_ylim{cluster_i});
    figure('Renderer', 'painters', 'Position', [100 100 300 500]);
    main_sdf_fig_out.draw
end

%% Get reference of neurons relative to analysis table
input_neurons = [];
input_neurons = find(mcc_analysis_table.glm_trial == 1);

for cluster_merge_i = 1:length(cluster_merge_idx)
    cluster_neuron_id_table{cluster_merge_i} = input_neurons(cluster_neuron_id{cluster_merge_i});
end
