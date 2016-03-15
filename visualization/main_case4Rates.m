

clear all;
clc;
close all;

%%

dirpath = fullfile('Data','Figure4');
name = 'Results_Results_VaryTransition_case3_';
binIndex = [5 10 15 20 25 30];
numBinIndex = length(binIndex);
Aindex = 1:7;
numIndex = length(Aindex);
numTrials = 5;

meanRate = zeros(numBinIndex,numIndex);
stdRate = zeros(numBinIndex,numIndex);
meanBinRate = zeros(numBinIndex,numIndex);
stdBinRate = zeros(numBinIndex,numIndex);
for q = 1:7
    for p = 1:numBinIndex
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

results.meanRate = meanRate;
results.stdRate = stdRate;
results.meanBinRate = meanBinRate;
results.stdBinRate = stdBinRate;
save('case4_rate.mat','results');

%%


load('case4_rate.mat');
meanRate = results.meanRate;
meanRate = [0 0.36 0.6, 1.2, 1.8, 2.4, 3.6];
stdRate = results.stdRate;
meanBinRate = results.meanBinRate;
stdBinRate = results.stdBinRate;

binIndex = [5 10 15 20 30];
plotIndex = [1 2 3 4 6];
Aindex = 1:7;
R = meanRate(1,:);

figure; hold on; box on;
legendname = '';
k = 1;
for i = Aindex
    errorbar(binIndex,meanBinRate(plotIndex,i),stdRate(plotIndex,i),'linewidth',2);
    legendname = [legendname '''R=' num2str(R(k),'%10.1f') ''','];
    k = k + 1;
end
eval(['h = legend(' legendname(1:end-1) ');']);
set(gca,'fontsize',28,'linewidth',2);
set(h,'fontsize',28,'location','northwest','box','off');
xlabel('Bin Size (steps)','fontsize',28);
ylabel('Transition Rate (per track)','fontsize',28);
axis([0 34 0 1]);
set(gca,'ytick',[0 .4 .8]);
print('-depsc','case4_binrate');

%% 

% transition rate vs bin size
% optimal bin size vs transition rate
% max percentage correct vs transition rate
% 

load('case4statistics.mat');
meanPercent = results2.meanPercent;
stdPercent = results2.stdPercent;
meanlogL = results2.meanLogL;
stdlogL = results2.stdLogL;

binIndex = [5 10 15 20 30];
plotIndex = [1 2 3 4 6];
Aindex = 1:7;
R = meanRate(1,:);


figure; hold on; box on;
for i = 1:7
    errorbar(binIndex,meanlogL(plotIndex,i),stdlogL(plotIndex,i),'linewidth',2)
end
set(gca,'fontsize',28,'linewidth',2);
xlabel('Bin Size (steps)','fontsize',28);
ylabel('Log-Likelihood','fontsize',28);
axis([0  34 1.32 1.65]);
set(gca,'ytick',[1.4 1.5 1.6]);
print('-depsc','case4_likelihood');


figure; hold on; box on;
for i = 1:7
    plot(binIndex,meanPercent(plotIndex,i),'linewidth',2)
end
set(gca,'fontsize',28,'linewidth',2);
xlabel('Bin Size (steps)','fontsize',28);
ylabel('Percent Correct','fontsize',28);
axis([0  34 0 1]);
print('-depsc','case4_percentcorrect');


[maxPercent,index] = max(meanPercent); 
error = [];
for i = 1:7
    error(i) = stdPercent(index(i),i);
end

figure; hold on; box on;
errorbar(R,maxPercent,error,'linewidth',2);
set(gca,'fontsize',28,'linewidth',2);
xlabel('Transition Rate (per track)','fontsize',28);
ylabel('Fraction Correct','fontsize',28);
axis([-.1  3.8 0.8 1]);
print('-depsc','case4_percentcorrect');


binIndex = [5 10 15 20 25 30];
figure; hold on; box on;
h = plot(R,binIndex(index),'linewidth',1);
scatter(R,binIndex(index),150,'markeredgecolor',h.Color,'markerfacecolor',h.Color,'linewidth',1.5);
set(gca,'fontsize',28,'linewidth',2);
xlabel('Transition Rate (per track)','fontsize',28);
ylabel('Optimal Bin (steps)','fontsize',28);
axis([-.1  3.8 14 31]);
print('-depsc','case4_optimalbin');

