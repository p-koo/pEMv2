

clear all;
clc;
close all;

%%

name = 'Results_case2';
binIndex = [5 10 15 20 25 30];
numBin = length(binIndex);

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



%%


clear all;
clc;
close all;

%%

name = 'Results_case2';
Nindex = [15 30 60 120 240];
numIndex = length(Nindex);

meanRate = zeros(numIndex,1);
stdRate = zeros(numIndex,1);
meanBinRate = zeros(numIndex,1);
stdBinRate = zeros(numIndex,1);
for p = 1:numIndex
    
    trialRate = zeros(numTrials,1);
    trialBinRate = zeros(numTrials,1);
    for z = 1:numTrials
        filename = [name '_' num2str(p) '_' num2str(z) '.mat'];
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









