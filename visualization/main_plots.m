clear all;
clc; 
close all;

%% load data 

% load pEMv2 results
[filename,dirpath] = uigetfile('*.mat','Select results file');
data = load(fullfile(dirpath,filename));
results = data.results;
splitX = results.X;
trackInfo = results.trackInfo;
splitLength = trackInfo.splitLength;
numFeatures = trackInfo.numFeatures;
vacf_exp = trackInfo.vacf_exp;
BIC = results.BIC;
state = results.state;

% load simulation data (i.e. ground truth)
[filename,dirpath] = uigetfile('*.mat','Select simulation file');
data = load(fullfile(dirpath,filename));
X = data.X;
markovStateSeq = data.markovStateSeq;
simParams = data.simParams;
true_numStates = simParams.numStates;
true_state = simParams.state;
true_A = simParams.A;


%% Print summary of results

disp('-------------------------------------------------------');
disp(['Summary of Results']);
disp(['True number of states: ' num2str(simParams.numStates) ]);
disp(['Est number of states: ' num2str(results.optimalSize) ]);
if results.optimalSize == simParams.numStates
    disp('Correct number of states found');
else
    disp('Wrong number of states found...');
    disp('WARNING: remaining analysis assumes correct number of states');
end
    
numStates = simParams.numStates;
est_vacf = state(simParams.numStates).vacf;
posteriorProb = state(simParams.numStates).posteriorProb;

%% Fraction correct (Ground truth state matches estimated state)

% true state sequence (mode within a bin)
true_stateIndex = zeros(length(splitX),1);
k = 1;
for i = 1:length(markovStateSeq)
    N = length(markovStateSeq{i});
    for j = 1:floor(N/splitLength)
        range = j*splitLength-splitLength+1:j*splitLength;
        true_stateIndex(k) = mode(markovStateSeq{i}(range));
        k = k + 1;
    end
end

% estimated state sequence
[MAX,est_stateIndex] = max(posteriorProb,[],2);
    
% figure out maximum matches by permuting state indices
permIndex = perms(1:numStates);
match = zeros(size(permIndex,1),1);
for j = 1:size(permIndex,1)
    tmpIndex = [];
    for i = 1:numStates
        index = find(est_stateIndex == permIndex(j,i));
        tmpIndex(index,1) = i;
    end
    match(j) = sum(true_stateIndex == tmpIndex)/length(true_stateIndex);
end
[MAX,index] = max(match);
fractionCorrect = MAX;
 
tmpIndex = zeros(length(true_stateIndex),1);
for i = 1:numStates
    tmpIndex(est_stateIndex == permIndex(index,i),1) = i;
end
est_stateIndex = tmpIndex;
match_Index = permIndex(index,:);

disp(['Fraction correct: ' num2str(fractionCorrect)]);


%% Compare ground truth and estimated MSD for each state

numLags = size(splitX{1},1)-1;
[true_stateMSD,true_stateMSDerror,true_stateN] = AverageMSD(splitX,true_stateIndex,numLags);
[est_stateMSD,est_stateMSDerror,est_stateN] = AverageMSD(splitX,est_stateIndex,numLags);

% plot MSD
figure; hold on; box on;
colorSet = hsv(numStates);
timeLags = 1:numLags;
legendname = 'h = legend(';
for i = 1:numStates
    stderror = est_stateMSDerror(i,:)./sqrt(est_stateN(i));
    errorbar(timeLags,est_stateMSD(i,:),stderror,'color',colorSet(i,:),'linewidth',1.5);
    legendname = [legendname '''State ' num2str(i) ''','];
end
for i = 1:numStates
    stderror = true_stateMSDerror(i,:)./sqrt(true_stateN(i));
    errorbar(timeLags,true_stateMSD(i,:),stderror,'--','color',colorSet(i,:),'linewidth',1.5);
end
set(gca,'fontsize',20,'linewidth',2);
eval([legendname(1:end-1) ');']);
set(h,'box','off','location','northwest','fontsize',20);
xlabel('Time lags (steps)','fontsize',20);
ylabel('MSD (\mum^2)','fontsize',20);


%% Compare ground truth and estimated covariance structure

[true_stateVACF,true_stateVACFerror,true_stateN] = AverageVACF(splitX,true_stateIndex);
[est_stateVACF,est_stateVACFerror,est_stateN] = AverageVACF(splitX,est_stateIndex);

% plot VACF
figure; hold on; box on;
colorSet = hsv(numStates);
timeLags = 1:numLags;
legendname = 'h = legend(';
for i = 1:numStates
    stderror = est_stateVACFerror(i,:)./sqrt(est_stateN(i));
    errorbar(timeLags,est_stateVACF(i,:),stderror,'color',colorSet(i,:),'linewidth',1.5);
    legendname = [legendname '''State ' num2str(i) ''','];
end
for i = 1:numStates
    stderror = true_stateVACFerror(i,:)./sqrt(true_stateN(i));
    errorbar(timeLags,true_stateVACF(i,:),stderror,'--','color',colorSet(i,:),'linewidth',1.5);
end
set(gca,'fontsize',20,'linewidth',2);
eval([legendname(1:end-1) ');']);
set(h,'box','off','location','northeast','fontsize',20);
xlabel('Time lags (steps)','fontsize',20);
ylabel('Covariance (\mum^2)','fontsize',20);


%% Compare ground truth and estimated transition matrix 

space = .6;
width = 30;
height = 50;

colorSet = 'rgbcm';

figure; hold on; box on;
k = 1;
for i = 1:height
    for j = 1:width
        x = X{k}(:,1) + i*space;
        y = X{k}(:,2) + j*space;
        plot(x,y,colorSet(trueState(k)));
        k = k + 1;
    end
end
axis equal;
axis([-1 height*space+1 -1.5 width*space+1]);
set(gca,'xtick',[],'ytick',[]);
rectangle('position',[-.5 -1 5 .3],'facecolor','k');
print('-depsc','case1_tracks_true');

figure; hold on; box on;
k = 1;
for i = 1:height
    for j = 1:width
        x = X{k}(:,1) + i*space;
        y = X{k}(:,2) + j*space;
        plot(x,y,colorSet(classifyState(k)));
        k = k + 1;
    end
end
axis equal;
axis([-1 height*space+1 -1.5 width*space+1]);
set(gca,'xtick',[],'ytick',[]);
print('-depsc','case1_tracks_classify');


%% plot each track (ground truth - left: estimated - right)



%% 





        