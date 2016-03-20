function [stateMSD,stateMSDerror,stateN] = AverageMSD(X,stateIndex,numLags)

numStates = max(stateIndex);
stateMSD = zeros(numStates,numLags);
stateMSDerror = zeros(numStates,numLags);
for i = 1:numStates
    index = find(stateIndex == i);
   
    stateN(i) = length(index);
    msd = zeros(stateN(i),numLags);
    for j = 1:stateN(i)
        msd(j,:) = CalculateMSD(X{index(j)},numLags);
    end
    stateMSD(i,:) = mean(msd);
    stateMSDerror(i,:) = std(msd);
end


