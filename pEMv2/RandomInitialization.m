function [vacf0,P0] = RandomInitialization(numStates,vacf_exp,method)
%--------------------------------------------------------------------------
% This function generates the model parameters, i.e. covariance matrix elements
% and population fractions, either via k-mean clustering or random
% initialization
%
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------


numFeatures = size(vacf_exp,2);
vacf_base = mean(vacf_exp,3);

switch method
    case 1
        % k-means clustering between the vacf of all particle tracks
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
        % random initialization (similar to pEM)
        
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
        vacf0 = zeros(numStates,numFeatures); 
        for i = 1:numStates
            index = find(cumprob > Dpoints(i),1,'first');
            C1 = edges(index);
            [vals, sort_index] = sort((vacf_base(:,1)-C1).*2, 'ascend');
            C2 = mean(vacf_base(sort_index(1:100),2),1);
            vacf0(i,1:2) = [C1 C2];
        end
        
    otherwise
end




