function [vacf0,P0] = RandomInitialization(numStates,vacf_exp,method)

numFeatures = size(vacf_exp,2);
vacf_base = mean(vacf_exp,3);

switch method
    case 1
        stateIndex = kmeans(vacf_base,numStates);
        
        vacf0 = zeros(numStates,numFeatures); P0 = zeros(1,numStates);
        for i = 1:numStates
            vacf0(i,:) = mean(vacf_base(stateIndex == i,:));
            P0(i) = sum(stateIndex == i);
        end
        P0 = P0/sum(P0);

        [vals,index] = sort(vacf0(:,1),'ascend');
        vacf0 = vacf0(index,:);
        P0 = P0(index);
    case 2
        
        P0 = zeros(1,numStates);
        for i = 1:numStates
            P0(i) = rand;
        end
        P0 = P0/sum(P0);

        bin = 1000;
        MIN = min(vacf_base(:,1));
        MAX = max(vacf_base(:,1));
        edges = MIN:(MAX-MIN)/bin:MAX;
        DN = histc(vacf_base(:,1),edges);
        cumprob = cumsum(DN)/sum(DN);

        % mean location of each random population fraction along CDF 
        Dpoints = cumsum([0 P0(1:end-1)] + diff([0 P0])/2);

        % find which diffusivity corresponds to each random population
        % fraction
        C1 = zeros(numStates,1);
        for i = 1:numStates
            C1(i) = edges(find(cumprob > Dpoints(i),1,'first'));
        end
        C2 = mean(vacf_base(:,2:end));
        
        vacf0 = zeros(numStates,numFeatures); 
        for i = 1:numStates
            vacf0(i,:) = [C1(i) C2];
        end
        
    otherwise
end

% CDF of diffusivities to generate a good random initialization



