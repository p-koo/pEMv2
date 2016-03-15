function [X,markovStateSeq] = SimulateDiffusionTransitions(simParams,numTracks,N,dt,numSubSteps)

% microstep times
microtime = dt/numSubSteps;
state = simParams.state;
A = simParams.A;
Pindex = simParams.Pindex;
locNoise = [state(:).locNoise]';
% A = 1;
% initial Markov state of each protein trajectory
markovStateSeq = SimulateMarkovState(numTracks,N,Pindex,A);

% simulate diffusion with transitions + static and dynamic localization noise
X = cell(numTracks,1);
for z = 1:numTracks    
%     disp(num2str(z));
    markovState = markovStateSeq{z};
    
    % simulate positions
    allTruePositions = zeros(N*numSubSteps,2);    
    
    lastState = markovState(1);
    lastPosition = [0 0];
    confineLocation = [];

    k = 1;
    while k < N
        lastIndex = find(markovState(k:end) ~= lastState,1,'first');
        if isempty(lastIndex)
            lastIndex = N;
        end
        range =  (k*numSubSteps-numSubSteps+1:(k+lastIndex-1)*numSubSteps);
        numSteps = length(range);

        D = state(lastState).D;
        alpha = state(lastState).A;
        v = state(lastState).v;
        L = state(lastState).L + randn(1)*state(lastState).Lnoise;        
        switch state(lastState).mode
            
            case 1 % normal diffusion
                deltaX = randn(numSteps,2).*sqrt(2*D*microtime);
                subPositions = cumsum([lastPosition; deltaX]);
            
            case 2 % confined diffusion    
                confineLocation = [confineLocation; lastPosition(1,1)-L lastPosition(1,1)+L lastPosition(1,2)-L lastPosition(1,2)+L];
                
                subPositions = zeros(numSteps+1,2);     
                subPositions(1,:) = lastPosition;
                for j = 2:numSteps + 1

                    % propose new position
                    proposedPosition = subPositions(j-1,:) + sqrt(2*D*microtime)*randn(1,2);

                    % impose reflecting boundary conditions
                    for i = 1:2
                        if proposedPosition(i) > confineLocation(end,i*2)
                            subPositions(j,i) = 2*confineLocation(end,i*2) - proposedPosition(i); 
                        elseif proposedPosition(i) < confineLocation(end,i*2-2+1)
                            subPositions(j,i) = 2*confineLocation(end,i*2-2+1) - proposedPosition(i); 
                        elseif proposedPosition(i) >= confineLocation(end,i*2-2+1) & proposedPosition(i) <= confineLocation(end,i*2)
                            subPositions(j,i) = proposedPosition(i);
                        end
                    end
                end
                
            case 3 % driven diffusion
                deltaX = randn(numSteps,2).*sqrt(2*D*microtime); 
                deltaX(:,1) = deltaX(:,1) + v*microtime;
                subPositions = cumsum([lastPosition; deltaX]);

            case 4 % fractional Brownian motion                
                model = 2*D;
                C = zeros(numSteps,numSteps);
                for n = 1:numSteps
                    for m = 1:numSteps
                        C(n,m) = model/2*(microtime)^alpha*(abs(n-m+1)^alpha - 2*abs(n-m)^alpha + abs(n-m-1)^alpha);
                    end
                end
                cholC = chol(C,'lower');
                
                deltaX = cholC*randn(numSteps,2);
                subPositions = cumsum([lastPosition; deltaX]);
        
            otherwise
        end
        
        % save true positions and state
        allTruePositions(range,:) = subPositions(2:end,:);
        
        if lastIndex == N
            k = N;
        else            
            lastState = markovState(k+lastIndex-1);
            k = k + lastIndex - 1;
            lastPosition = subPositions(end,:);
        end
    end

    % find true positions every dt
    truePositions = allTruePositions(1:numSubSteps:N*numSubSteps,:);

    % average over true positions to get the motion blur 
    motionPositions = zeros(N,2);
    for i = 1:N
        range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
        motionPositions(i,:) = mean(allTruePositions(range,:));
    end

    % simulate noise 
    positionNoise = randn(N,2).*repmat(locNoise(markovState),1,2);

    % get the observed displacements
    X{z} = motionPositions + positionNoise;
    
%     figure(1); clf;
%     plot(X{z}(:,1),X{z}(:,2));
%     pause;
end
    





%%

    
    
    
% function [X,markovStateSeq] = SimulateDiffusionTransitions(simParams,numTracks,N,dt,numSubSteps)
% 
% % microstep times
% microtime = dt/numSubSteps;
% state = simParams.state;
% A = simParams.A;
% Pindex = simParams.Pindex;
% locNoise = [state(:).locNoise]';
% 
% % calculate covariance matrix if fBM
% cholC = cell(length(state),1);
% for i = 1:length(state)
%     model = 2*state(i).D;
%     alpha = state(i).A;
%     if alpha ~= 1
%         C = zeros(numSubSteps,numSubSteps);
%         for n = 1:numSubSteps
%             for m = 1:numSubSteps
%                 C(n,m) = model/2*(microtime)^alpha*(abs(n-m+1)^alpha - 2*abs(n-m)^alpha + abs(n-m-1)^alpha);
%             end
%         end
%         cholC{i} = chol(C,'lower');
%     end
% end
% 
% % initial Markov state of each protein trajectory
% Nindex = round(numTracks*Pindex);
% initialState = [];
% for i = 1:simParams.numStates
%     initialState = [initialState ones(1,Nindex(i))*i];
% end
% if length(initialState) < numTracks
%     initialState = [ones(1,numTracks-length(initialState)) initialState];
% end
% initialMarkovState = initialState;
% 
% % start location of each protein trajectory
% startLocation = zeros(numTracks,1);
% 
% % map transition probabilities onto uniform distribution
% transitionCum = cumsum(A,2);
% 
% % simulate diffusion with transitions + static and dynamic localization
% % noise
% X = cell(numTracks,1); markovStateSeq = cell(numTracks,1); 
% for z = 1:numTracks
% 
%     % simulate transition probabilities for each protein trajectory
%     markovStateChange = rand(N,1);
% 
%     % underlying markov markovState sequences
%     markovState = zeros(N,1);
%     markovState(1) = initialMarkovState(z);
%     for t = 2:N
%         markovState(t) = find(markovStateChange(t) < transitionCum(markovState(t-1),:),1);
%     end
%     
%     % set up initial confinement
%     switch state(markovState(1)).mode
% 
%         case 1 % normal diffusion
%             confineLocation = [];
% 
%         case 2 % confined diffusion                
%             L = state(markovState(1)).L + randn(1)*state(markovState(1)).Lnoise;
%             confineLocation = [-L L -L L];
%         otherwise
% 
%     end
% 
%     % simulate positions
%     lastState = markovState(1);
%     allTruePositions = zeros(N*numSubSteps,2);
%     allmarkovStateSeq = zeros(N*numSubSteps,2);
%     allTruePositions(1,:) = startLocation(z,:);
%     for t = 1:N
%         range = 1 + (t*numSubSteps-numSubSteps+1:t*numSubSteps);
% 
%         % free diffusion
%         D = state(markovState(t)).D;
%         switch state(markovState(t)).mode
%             
%             case 1 % normal diffusion
%                 deltaX = randn(numSubSteps,2).*sqrt(2*repmat(D,numSubSteps,2)*microtime);
%                 subPositions = cumsum([allTruePositions(range(1)-1,:); deltaX]);
%             
%             case 2 % confined diffusion                
%                 L = state(markovState(t)).L + randn(1)*state(markovState(t)).Lnoise;
%                 
%                 lastPosition = allTruePositions(range(1)-1,:);
% 
%                 if markovState(t) ~= lastState
%                     confineLocation = [confineLocation; lastPosition(1,1)-L lastPosition(1,1)+L lastPosition(1,2)-L lastPosition(1,2)+L];
%                 end
% 
%                 subPositions = zeros(numSubSteps+1,2);
%                 subPositions(1,:) = lastPosition;
%                 for j = 2:numSubSteps + 1
% 
%                     % propose new position
%                     proposedPosition = subPositions(j-1,:) + sqrt(2*D*microtime)*randn(1,2);
% 
%                     % impose reflecting boundary conditions
%                     for i = 1:2
%                         if proposedPosition(i) > confineLocation(end,i*2)
%                             subPositions(j,i) = 2*confineLocation(end,i*2) - proposedPosition(i); 
%                         elseif proposedPosition(i) < confineLocation(end,i*2-2+1)
%                             subPositions(j,i) = 2*confineLocation(end,i*2-2+1) - proposedPosition(i); 
%                         elseif proposedPosition(i) >= confineLocation(end,i*2-2+1) & proposedPosition(i) <= confineLocation(end,i*2)
%                             subPositions(j,i) = proposedPosition(i);
%                         end
%                     end
%                 end
%                 
%             case 3 % driven diffusion
%                 deltaX = randn(numSubSteps,2).*sqrt(2*repmat(D,numSubSteps,2)*microtime) + state(markovState(t)).v*microtime;
%                 subPositions = cumsum([allTruePositions(range(1)-1,:); deltaX]);
% 
%             case 4 % fractional Brownian motion
%                 deltaX = cholC{markovState(t)}*randn(numSubSteps,2);
%                 subPositions = cumsum([allTruePositions(range(1)-1,:); deltaX]);
%                 
%             otherwise
%         end
%         
%         % set new last markov state
%         lastState = markovState(t);
%         
%         % save true positions and state
%         allTruePositions(range,:) = subPositions(2:end,:);
%         allmarkovStateSeq(range) = ones(length(range),1)*lastState;
%     end
% 
%     % find true positions every dt
%     truePositions = allTruePositions(1:numSubSteps:N*numSubSteps,:);
% 
%     % average over true positions to get the motion blur 
%     motionPositions = zeros(N,2);
%     for i = 1:N
%         range = i*numSubSteps-numSubSteps+1:i*numSubSteps;
%         motionPositions(i,:) = mean(allTruePositions(range,:));
%     end
% 
%     % simulate noise 
%     positionNoise = randn(N,2).*repmat(locNoise(markovState),1,2);
% 
%     % get the observed displacements
%     X{z} = motionPositions + positionNoise;
%     markovStateSeq{z} = markovState;
% end
% 
% 
% 
