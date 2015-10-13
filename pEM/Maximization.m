function [vacf,P] = Maximization(gamma,vacf_exp)

% get track population parameters
[numTracks,numStates,dim] = size(gamma);
numFeatures = size(vacf_exp,2);

% population fraction maximization
Nk = zeros(dim,numStates);
for i = 1:dim
    Nk(i,:) = sum(gamma(:,:,i));
end
meanNk = mean(Nk,1);
P = meanNk/numTracks;

% diffusion coefficient maximization
vacf = zeros(numStates,numFeatures);
for j = 1:dim
    for k = 1:numFeatures
        vacf(:,k) = vacf(:,k) + (sum(gamma(:,:,j).*(vacf_exp(:,k,j)*ones(1,numStates)))./Nk(j,:))';
    end
end
vacf = vacf/dim;

