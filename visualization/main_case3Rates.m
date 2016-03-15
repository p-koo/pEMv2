

clear all;
clc;
close all;

%%

dirpath = fullfile('Data','Figure3');
name = 'Results_VaryBin_case2_5_';
binIndex = [5 10 15 20 25 30];
numBin = length(binIndex);
numTrials = 5;

meanRate = zeros(numBin,1);
stdRate = zeros(numBin,1);
meanBinRate = zeros(numBin,1);
stdBinRate = zeros(numBin,1);
for p = 1:numBin
    binSize = binIndex(p);

    trialRate = zeros(numTrials,1);
    trialBinRate = zeros(numTrials,1);
    for z = 1:numTrials
        filename = [name num2str(z) '_' num2str(binSize) '.mat'];
        load(fullfile(dirpath,filename));

        stateSeq = results.stateSeq;
        numTracks = length(stateSeq);
        numBin = length(results.binIndex);

        numTransitions = 0;
        for i = 1:numTracks
            seq = stateSeq{i};
            while length(seq) > 0
                state = seq(1);
                index = find(seq ~= state,1,'first');
                if ~isempty(index)
                    numTransitions = numTransitions + 1;
                    seq(1:index-1) = [];
                else
                    break;
                end
            end
        end
        trialRate(z) = numTransitions/numTracks;
        trialBinRate(z) = numTransitions/numBin;
    end

    meanRate(p) = mean(trialRate);
    stdRate(p) = std(trialRate);
    meanBinRate(p) = mean(trialBinRate);
    stdBinRate(p) = std(trialBinRate);
end


results.meanRate = meanRate;
results.stdRate = stdRate;
results.meanBinRate = meanBinRate;
results.stdBinRate = stdBinRate;
save('case3_rate.mat','results');


%%

load('case3_rate.mat');
meanRate = results.meanRate;
stdRate = results.stdRate;
meanBinRate = results.meanBinRate;
stdBinRate = results.stdBinRate;

binIndex = [5 10 15 20 30];
plotIndex = [1 2 3 4 6];

figure; hold on; box on;
errorbar(binIndex,meanBinRate(plotIndex),stdRate(plotIndex),'linewidth',2);
set(gca,'fontsize',28,'linewidth',2);
xlabel('Bin Size (steps)','fontsize',28);
ylabel('Transition Rate (per track)','fontsize',28);
set(gca,'ytick',[0 .1 .2 .3]);
axis([4.5 30.5  0 .33]);
print('-depsc','case3_rate');









