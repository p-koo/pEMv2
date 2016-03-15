clear all;
clc;
close all;

%%

savepath = fullfile('Results','VaryN');

% simulation parameters
simSet = 2;
uniformA = 0;
dt = .055;
numRepeat = 5;

% analysis parameters
minStates = 2;
maxStates = 4;
correctDrift = 1;           % drift velocity correction
splitLength = 10;           % bin size of protein trajectories
numFeatures = 9;            % number of features (must be less than bin size)
numReinitialize = 25;       % number of reinitialization trials
numPerturbation = 0;        % number of perturbations trials
toleranceEM = 1e-7;         % convergence condition for EM
maxiterEM = 100000;         % maximum number of iterations for EM
verboseEM = 0;              % display progress of algorithm
maxiterHMM = 100;           % maximum number of iterations for EM
toleranceHMM = 1e-5;        % convergence condition for EM
verboseHMM = 0;             % display progress of algorithm

%% BIC score

Nindex = [10 30 60 120 240 480];
trackIndex = [4800 1600 800 400 200 100];
for j = 1:length(Nindex)
    N = Nindex(j);
    numTracks = trackIndex(j);
    disp(['N = ' num2str(N)]);    
    
    load(fullfile('data',['VaryN_case2_' num2str(j) '.mat']));
    load(fullfile(savepath,['Case2_N=' num2str(N)]));

    
    % track specifications
    trackInfo.splitLength = splitLength;
    trackInfo.numFeatures = numFeatures;
    trackInfo.drift = correctDrift;
    trackInfo.rotate = correctDrift;
    trackInfo.dt = dt;             
    trackInfo.R = 1/6;             

    % EM parameters
    emparams.numReinitialize = numReinitialize;             
    emparams.numPerturbation = numPerturbation;    
    emparams.tolerance = toleranceEM;             
    emparams.maxiter = maxiterEM;               
    emparams.verbose = verboseEM;

    % HMM parameters
    hmmparams.maxiter = maxiterHMM;
    hmmparams.tolerance = toleranceHMM;
    hmmparams.verbose = verboseHMM;

    
    for z = 1:numRepeat
        X = tracks(z).X;
        stateSeq = tracks(z).stateSeq;

        % split tracks into equally spaced bin sizes
        [splitX,splitState,splitIndex,trackInfo] = SplitTracks(X,stateSeq,trackInfo);
        logL = zeros(1,maxStates);
        for numStates = minStates:maxStates
            
            gammank = trial(z).results(numStates).hmm.gammank;
            mu = trial(z).results(numStates).hmm.mu;
            sigma = trial(z).results(numStates).hmm.sigma;
            vacf = zeros(numStates,size(sigma,1));
            for i = 1:numStates
                vacf(i,:) = sigma(1,:,i);
            end
            P = sum(gammank);
            P = P/sum(P);
            [gamma,L] = Expectation(splitX,vacf,mu,P,trackInfo);
             
%             p = trial(z).results(numStates).hmm.p;
%             a = trial(z).results(numStates).hmm.a;
%             b = trial(z).results(numStates).hmm.b;
%             [est,L,gammank] = ExpectationHMM(p,a,b,splitIndex);
%             [p,a,b,mu,sigma] = MaximizationHMM2(splitX,gammank,est,numStates,trackInfo);
%             [est,L,gammank] = ExpectationHMM(p,a,b,splitIndex);
            logL(numStates) = -L;
%             logL2(z,numStates) = -trial(z).results(numStates).hmm.logL;
        end
        
        BIC = -logL - (1:maxStates)*log(length(splitX))
        [MAX index] = max(BIC(2:end));
        disp(num2str(index+1));
    end
    
end


%% percent match

percentMatch = zeros(length(Nindex),numRepeat);
for j = 1:length(Nindex)
    N = Nindex(j);
    numTracks = trackIndex(j);
    disp(['N = ' num2str(N)]);        
    load(fullfile(savepath,['Case2_N=' num2str(N)]));

    for z = 1:numRepeat        
        percentMatch(j,z) = trial(z).compare.percentMatch;
    end
end

figure; errorbar(Nindex,mean(percentMatch,2),std(percentMatch,[],2))

%% log-loss

logLoss = zeros(length(Nindex),numRepeat);
for j = 1:length(Nindex)
    N = Nindex(j);
    numTracks = trackIndex(j);
    disp(['N = ' num2str(N)]);        
    load(fullfile(savepath,['Case2_N=' num2str(N)]));

    for z = 1:numRepeat        
        logLoss(j,z) = mean(trial(z).compare.logLoss);
    end
end
figure; errorbar(Nindex,mean(logLoss,2),std(logLoss,[],2))


%% transition rate

trueRate = zeros(length(Nindex),numRepeat);
stateRate = zeros(length(Nindex),numRepeat);
for j = 1:length(Nindex)
    N = Nindex(j);
    numTracks = trackIndex(j);
    disp(['N = ' num2str(N)]);        
    load(fullfile(savepath,['Case2_N=' num2str(N)]));

    for z = 1:numRepeat        
        trueRate(j,z) = mean(trial(z).compare.trueRate);
        stateRate(j,z) = mean(trial(z).compare.stateRate);
    end
end

figure; hold on;
errorbar(Nindex,mean(trueRate,2),std(trueRate,[],2),'--k')
errorbar(Nindex,mean(stateRate,2),std(stateRate,[],2),'r')

%% msd comparison

numLags = 25;
trueMSD = cell(1,3);
stateMSD = cell(1,3);
errortrueMSD = cell(1,3);
errorstateMSD = cell(1,3);
for i = 1:3
    trueMSD{i} = zeros(length(Nindex),numLags);
    stateMSD{i} = zeros(length(Nindex),numLags);
    errortrueMSD{i} = zeros(length(Nindex),numLags);
    errorstateMSD{i} = zeros(length(Nindex),numLags);

    for j = 1:length(Nindex)
        N = Nindex(j);
        numTracks = trackIndex(j);
        load(fullfile(savepath,['Case2_N=' num2str(N)]));

        tmpMSD1 = zeros(numRepeat,numLags);
        tmpMSD2 = zeros(numRepeat,numLags);
        for z = 1:numRepeat        
            tmpMSD1(z,:) = trial(z).compare.truemsd(i,:);
            tmpMSD2(z,:) = trial(z).compare.statemsd(i,:);
        end
        trueMSD{i}(j,:) = mean(tmpMSD1);
        stateMSD{i}(j,:) = mean(tmpMSD2);
        errortrueMSD{i}(j,:) = std(tmpMSD1);
        errorstateMSD{i}(j,:) = std(tmpMSD2);
    end
end

color = 'bgrmc';
close all;
for j = 1:length(Nindex)
    figure; hold on;
    for i = 1:3
        errorbar((1:numLags)*dt,trueMSD{i}(j,:),errortrueMSD{i}(j,:),['--' color(i)]);
        errorbar((1:numLags)*dt,stateMSD{i}(j,:),errorstateMSD{i}(j,:),color(i));
    end
end








