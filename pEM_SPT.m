function  results = pEM_SPT(X,trackInfo,params)

minStates = params.minStates;
maxStates = params.maxStates;
numFeatures = params.numFeatures;
dt = trackInfo.dt;
numData = length(X)*(length(X{1})-1);

% calculate the displacements for each particle track
deltaX = cell(trackInfo.numberOfTracks,1);
for i = 1:trackInfo.numberOfTracks
    deltaX{i} = diff(X{i});
end

% calculate relevant properties to enhance compuatational time
[trackInfo.vacf_exp,trackInfo.xbar_exp] = CovarianceProperties(deltaX,numFeatures);

% BIC Model Selection Loop
state = struct([]);
BIC = zeros(maxStates,1); 
for numStates = minStates:maxStates
    disp([num2str(numStates) ' state model']);

    % random initialization
    [vacf0,P0] = RandomInitialization(numStates,trackInfo.vacf_exp,1);
    if params.verbose == 1
        disp(['Initial Guess: ' num2str((vacf0(:,1) + 2*vacf0(:,2))'/2/dt)]);
    end
    
    % run rEM
    [baseVacf,baseP,posteriorProb,logLmax,remTrial] = rEM(deltaX,vacf0,P0,params,trackInfo);    
    if params.verbose == 1
        disp(['rEM: ' num2str((baseVacf(:,1) + 2*baseVacf(:,2))'/2/dt)]);
    end
    
    % run pEM
    [baseVacf,baseP,posteriorProb,logLmax,pemTrial] = pEM(deltaX,baseVacf,baseP,logLmax,params,trackInfo);
    if params.verbose == 1
        disp(['pEM: ' num2str((baseVacf(:,1) + 2*baseVacf(:,2))'/2/dt)]);
    end
    
    % calculate BIC
    nparams = numStates + numStates*trackInfo.numFeatures; 
    BIC(numStates) = logLmax/2 - nparams/2*log(numData);
    
    % store results
    state(numStates).numberOfStates = numStates;
    state(numStates).BIC = BIC(numStates);
    state(numStates).logL = logLmax;
    state(numStates).vacf = baseVacf;
    state(numStates).P = baseP;
    state(numStates).posteriorProb = posteriorProb;
    state(numStates).remTrial = remTrial;
    state(numStates).pemTrial = pemTrial;
    
    if BIC(numStates) < max(BIC(numStates))
        break;
    end
end

% store optimal results
[MAX,numStates] = max(BIC);
results.X = X;
results.params = params;
results.trackInfo = trackInfo;
results.state = state;
results.optimalSize = numStates;
results.optimalVacf = state(numStates).vacf;
results.optimalP = state(numStates).P;
results.optimalL = state(numStates).logL;
results.posteriorProb = state(numStates).posteriorProb;
results.BIC = BIC;
