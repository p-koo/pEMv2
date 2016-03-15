

clear all;
clc;
close all;

%%

name = 'Results_Results_VaryTransition_case3_';
binIndex = [5 10 15 20 25 30];
numBin = length(binIndex);
Aindex = 1:7;
numIndex = length(Aindex);

meanRate = zeros(numBin,numIndex);
stdRate = zeros(numBin,numIndex);
meanBinRate = zeros(numBin,numIndex);
stdBinRate = zeros(numBin,numIndex);
for q = 1:7
    for p = 1:numBin
        binSize = binIndex(p);
        
        trialRate = zeros(numTrials,1);
        trialBinRate = zeros(numTrials,1);
        for z = 1:numTrials
            filename = [name num2str(q) '_' num2str(z) '_' num2str(binSize) '.mat'];
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
        
        meanRate(p,q) = mean(trialRate);
        stdRate(p,q) = std(trialRate);
        meanBinRate(p,q) = mean(trialBinRate);
        stdBinRate(p,q) = std(trialBinRate);
    end
end



