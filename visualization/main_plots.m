clear all;
clc; 
close all;

%% load data 

% load simulation data (i.e. ground truth)
[filename,dirpath] = uigetfile('*.mat','Select simulation file');
data = load(fullfile(dirpath,filename));
X = data.X;
true_StateSeq = data.markovStateSeq;
simParams = data.simParams;
true_numStates = simParams.numStates;
true_state = simParams.state;
true_A = simParams.A;


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
numTracks = length(true_StateSeq);
est_vacf = state(simParams.numStates).vacf;
posteriorProb = state(simParams.numStates).posteriorProb;


analytics = struct;


%% Fraction correct (Ground truth state matches estimated state)

% true state sequence (mode within a bin)
true_stateIndex = zeros(length(splitX),1);
k = 1;
for i = 1:length(true_StateSeq)
    N = length(true_StateSeq{i});
    for j = 1:floor(N/splitLength)
        range = j*splitLength-splitLength+1:j*splitLength;
        true_stateIndex(k) = mode(true_StateSeq{i}(range));
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
 
% reconstruct indices based on highest match
tmpIndex = zeros(length(true_stateIndex),1);
for i = 1:numStates
    tmpIndex(est_stateIndex == permIndex(index,i),1) = i;
end
est_stateIndex = tmpIndex;
match_Index = permIndex(index,:);

disp(['Fraction correct: ' num2str(fractionCorrect)]);


analytics.fractionCorrect = fractionCorrect;
analytics.true_stateIndex = true_stateIndex;
analytics.est_stateIndex = est_stateIndex;

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


analytics.true_stateMSD = true_stateMSD;
analytics.est_stateMSD = est_stateMSD;
analytics.true_stateMSDerror = true_stateMSDerror;
analytics.est_stateMSDerror = est_stateMSDerror;


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


analytics.true_stateVACF = true_stateVACF;
analytics.est_stateVACF = est_stateVACF;
analytics.true_stateVACFerror = true_stateVACFerror;
analytics.est_stateVACFerror = est_stateVACFerror;
analytics.true_stateN = true_stateN;
analytics.est_stateN = est_stateN;


%% Compare ground truth and estimated transition matrix 


% transition probability (true empirical)
A = zeros(numStates,numStates);
for i = 1:numTracks
    seq = true_StateSeq{i};    
    for n = 1:length(seq)-1
        A(seq(n),seq(n+1)) = A(seq(n),seq(n+1)) + 1;
    end
end
norm = sum(A,2);
true_transitionMatrix = A./repmat(norm,1,numStates);
        
% transition probability (classified)
est_StateSeq = cell(numTracks,1);
for i = 1:numTracks
    index = find(splitIndex == i);
    tmp_states = [];
    for j = 1:length(index)
        tmp_states = [tmp_states; ones(splitLength,1)*est_stateIndex(index(j))];
    end
    est_StateSeq{i} = tmp_states;
end

A = zeros(numStates,numStates);
for i = 1:length(est_StateSeq)
    seq = est_StateSeq{i};    
    for n = 1:length(seq)-1
        A(seq(n),seq(n+1)) = A(seq(n),seq(n+1)) + 1;
    end
end
norm = sum(A,2);
est_transitionMatrix = A./repmat(norm,1,numStates);

disp('true transition matrix:');
disp(num2str(true_transitionMatrix));
disp('est transition matrix:');
disp(num2str(est_transitionMatrix));


analytics.true_transitionMatrix = true_transitionMatrix;
analytics.est_transitionMatrix = est_transitionMatrix;

analytics.est_StateSeq = est_StateSeq;
analytics.true_StateSeq = true_StateSeq;

%% plot each track (ground truth - left: estimated - right)

offset = .05;
colors = 'bgrcmyk';
for i = 1:numTracks

    % classified states
    estSeq = est_StateSeq{i};
    N = length(estSeq);
    
    % true states
    trueSeq = true_stateSeq{i};
    trueSeq = trueSeq(1:N);
    
    x = X{i}(1:N,:);
    minx = min(x);
    maxx = max(x);
    
    figure(1); clf; 
    subplot(1,2,1); box on; hold on; axis equal;
    for j = 1:N-1
        plot(x(j:j+1,1),x(j:j+1,2),'color',colors(trueSeq(j)));
    end
    axis([minx(1)-offset maxx(1)+offset minx(2)-offset maxx(2)+offset]);
    set(gca,'linewidth',1.5,'xticklabel',[],'yticklabel',[]);
    title('Ground truth','fontsize',20);
    
    subplot(1,2,2); box on; hold on; axis equal;
    for j = 1:N-1
        plot(x(j:j+1,1),x(j:j+1,2),'color',colors(estSeq(j)));
    end    
    axis([minx(1)-offset maxx(1)+offset minx(2)-offset maxx(2)+offset]);
    set(gca,'linewidth',1.5,'xticklabel',[],'yticklabel',[]);
    
    % error bar is 200 nm
    rectangle('position',[minx(1) minx(2) .2 .02],'facecolor','k','edgecolor','k');
    title('Classified','fontsize',20);
    title(['track ' num2str(i)]);
    
    pause;
    
end



%% 





        