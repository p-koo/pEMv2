function [X] = SimulateDiffusionConst(simParams,numTracks,N,dt,numSubSteps)

% microstep times
microtime = dt/numSubSteps;
state = simParams.state;
locNoise = [state(:).locNoise]';

% simulate diffusion with transitions + static and dynamic localization noise
numSteps = N*numSubSteps;
model = 2*simParams.state.D;
alpha = simParams.state.A;
C = zeros(numSteps,numSteps);
for n = 1:numSteps
    for m = 1:numSteps
        C(n,m) = model/2*(microtime)^alpha*(abs(n-m+1)^alpha - 2*abs(n-m)^alpha + abs(n-m-1)^alpha);
    end
end
cholC = chol(C,'lower');

X = cell(numTracks,1);
for z = 1:numTracks    
    
    % simulate positions
    allTruePositions = zeros(N*numSubSteps,2);    
    
    deltaX = cholC*randn(numSteps,2);
    allTruePositions = cumsum([0 0; deltaX]);

    % find true positions every dt
    truePositions = allTruePositions(1:numSubSteps:N*numSubSteps,:);

    % average over true positions to get the motion blur 
    motionPositions = zeros(N,2);
    for i = 1:N
        range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
        motionPositions(i,:) = mean(allTruePositions(range,:));
    end

    % simulate noise 
    positionNoise = randn(N,2).*repmat(locNoise,N,2);

    % get the observed displacements
    X{z} = motionPositions + positionNoise;
    
%     figure(1); clf;
%     plot(X{z}(:,1),X{z}(:,2));
%     pause;
end
    