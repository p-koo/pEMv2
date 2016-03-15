function [statemsd,stateRate] = StateMSD(numLags,X,splitIndex,estIndex,splitLength)

numTracks = length(X);
numStates = max(estIndex);

expSeq = cell(numTracks,1);
for i = 1:numTracks
    index = find(splitIndex == i);
    expSeq{i} = zeros(length(X{i}),1);
    for j = 1:length(index)
        range = j*splitLength-splitLength+2-j:j*splitLength+1-j;
        expSeq{i}(range) = ones(length(range),1)*estIndex(index(j));
    end
    index = find(expSeq{i} == 0,1,'first');
    expSeq{i}(index:end) = [];
end

[stateX,numTransitions2]= ParseStateSeq(X,expSeq,numStates);
statemsd = AverageMSD(stateX,numStates,numLags);
stateRate = sum(numTransitions2-1)/numTracks;


