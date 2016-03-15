
load(fullfile('data',['VaryFeatures_case1.mat']));
savepath = fullfile('Results','VaryFeatures');
% load 

simSet = 1;
dt = .055;
N = 120;
numTracks = 3000;
numRepeat = 5;

%%

numStates = 4;
z = 1;
X = tracks(z).X;
stateSeq = tracks(z).stateSeq;
results = trial(z).results;
trueIndex = trial(z).compare.trueIndex;
estIndex = trial(z).compare.estIndex;
percentMatch = trial(z).compare.percentMatch;
        
a = results(numStates).hmm.a;
sigma = results(numStates).hmm.sigma;
mu = results(numStates).hmm.mu;
gammank = results(numStates).hmm.gammank;
logL = results(numStates).hmm.logL;

sum(trueIndex == estIndex)/length(trueIndex)

% calculate ground truth MSD
numLags = 25;
[truthX,numTransitions] = ParseStateSeq(X,stateSeq,numStates);
truemsd = AverageMSD(truthX,numStates,numLags);
trueRate = sum(numTransitions-1)/length(stateSeq);
disp(['Transitions per track: ' num2str(trueRate)]);

% calculate state MSD
[statemsd,stateRate] = StateMSD(numLags,X,splitIndex,estIndex,splitLength);

% compare MSD
color = 'rbgmcy';
figure; hold on; box on;
for i = 1:numStates
    plot(statemsd(i,:),color(i));
    plot(truemsd(i,:),['--' color(i)]);
end
set(gca,'yscale','log','fontsize',20);
xlabel('Time lags','fontsize',20);
ylabel('MSD','fontsize',20);
% print('-depsc','2state_bin=3');


%%

logL = zeros(maxStates,1);
BIC = zeros(maxStates,1);
for numStates = minStates:maxStates
    nparams = numStates.*(numStates-1) + numStates.*trackInfo.trackLength.*(trackInfo.trackLength-1) + numStates;
    BIC(numStates) = results(numStates).hmm.logL - nparams/2*log(trackInfo.numberOfTracks);
    logL(numStates) = results(numStates).hmm.logL;
end


%%

gamma = results(numStates).hmm.gammank;
gamma = gamma + 1e-3;
gamma = gamma./repmat(sum(gamma,2),1,numStates);
gamma = gamma(:,sIndex);


% index = [1 2 4 3];
% gamma = gamma(:,index);

logLoss = [];
for i = 1:length(stateIndex)
    logLoss = [logLoss; -log(gamma(i,stateIndex(i)))];
end
mean(logLoss)

figure; hold on;
[f x] = ecdf(logLoss);
plot(x,f,'k');
plot(-log(1/numStates)*ones(1,2),[0 1],'--r');
set(gca,'xscale','log');



