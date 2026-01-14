figure
hold on
colors = [[0 0.4470 0.7410 .4];[0.8500 0.3250 0.0980 .4];[0.4660 0.6740 0.1880 .4]]
compInds = [1 5 9];
for jj = 1:3
subfig = compInds(jj);
    for ii = 1:length(onsets{subfig})
        line([onsets{subfig}(ii)-1000, onsets{subfig}(ii)+ durs{subfig}(ii)-1000 ],...
            [AUC{subfig}(ii) AUC{subfig}(ii)],'color',colors(jj,:),'linewidth',2)
    end
end
compInds = [2 6 10];
for jj = 1:3
subfig = compInds(jj);
    for ii = 1:length(onsets{subfig})
        line([onsets{subfig}(ii)-1000, onsets{subfig}(ii)+ durs{subfig}(ii)-1000 ],...
            [AUC{subfig}(ii) AUC{subfig}(ii)],'color',colors(jj,:),'linewidth',2)
    end
end
xlim([-100 2100])
ylabel(['AUROC'])
xlabel(['Time from saccade (ms)'])
legend({'SEF','dACC','vACC'},'Location','east')

figure
hold on
colors = [[0 0.4470 0.7410 .4];[0.8500 0.3250 0.0980 .4];[0.4660 0.6740 0.1880 .4]];
compInds = [1 5 9];
for jj = 1:3
subfig = compInds(jj);
    for ii = 1:length(onsets{subfig})
        line([0 durs{subfig}(ii) ],...
            [AUC{subfig}(ii) AUC{subfig}(ii)],'color',colors(jj,:),'linewidth',2)
    end
end
compInds = [2 6 10];
for jj = 1:3
subfig = compInds(jj);
    for ii = 1:length(onsets{subfig})
        line([0 durs{subfig}(ii)],...
            [AUC{subfig}(ii) AUC{subfig}(ii)],'color',colors(jj,:),'linewidth',2)
    end
end
xlim([-100 2100])
ylabel(['AUROC'])
xlabel(['Tmodulation duration (ms)'])
legend({'SEF','dACC','vACC'},'Location','east')