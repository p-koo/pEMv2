function [stateVACF,stateVACFerror,stateN] = AverageVACF(X,stateIndex)

[N,dim] = size(X{1});
numStates = max(stateIndex);
stateVACF = zeros(numStates,N-1);
stateVACFerror = zeros(numStates,N-1);
for i = 1:numStates
    index = find(stateIndex == i);
   
    stateN(i) = length(index);
    vacf = zeros(stateN(i),N-1);
    for j = 1:stateN(i)
        C = 0;
        for k = 1:dim
            x = diff(X{index(j)}(:,k));
            C = C + x*x';
        end
        C = C/dim;
        vacf(j,:) = CalculateVACF(C);
    end
    stateVACF(i,:) = mean(vacf);
    stateVACFerror(i,:) = std(vacf);
end

