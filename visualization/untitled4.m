
binSizeIndex = [5 10 15 20 30];
plotRange = [1 2 3 4 6];
dirpath = fullfile('Data','Figure3');
name = ['Results_VaryBin_case2_'];
numIndex = length(binSizeIndex);
numTrials = 5;

for p = 1:numIndex
    binSize = binSizeIndex(p);
    for i = trialIndex
        filename = [name num2str(5) '_' num2str(i) '_' num2str(binSizeIndex(p)) '.mat'];
        data = load(fullfile(dirpath,filename));        
        X = data.results.X;
        state = data.results.state;
        stateSeq = data.results.stateSeq;
        binSize = size(results(k).trial(i).binX{1},1)+1;
        [binStateIndex,stateIndex,transitionIndex] = BinStateSequence(stateSeq,binSize,1);
        gamma = mean(state(numStates).pem.gamma,3);
        [MAX,estState] = max(gamma,[],2);
        
        
        
        
    end
    k = k + 1;
end

% transition probability (true empirical)
numStates = 3;
A = zeros(numStates,numStates);
for i = 1:length(stateSeq)
    seq = stateSeq{i};    
    for n = 1:length(seq)-1
        A(seq(n),seq(n+1)) = A(seq(n),seq(n+1)) + 1;
    end
end
norm = sum(A,2);
prob = A./repmat(norm,1,numStates);
        
% transition probability (classified)
[MAX,estState] = max(gamma,[],2);
estStateIndex = cell(length(stateSeq),1);
for i = 1:length(stateSeq)
    index = find(stateIndex == i);
    states = [];
    for j = 1:length(index)
        states = [states; ones(binSize,1)*estState(index(j))];
    end
    estStateIndex{i} = states;
end

A = zeros(numStates,numStates);
for i = 1:length(estStateIndex)
    seq = estStateIndex{i};    
    for n = 1:length(seq)-1
        A(seq(n),seq(n+1)) = A(seq(n),seq(n+1)) + 1;
    end
end
norm = sum(A,2);
prob = A./repmat(norm,1,numStates);



        
        
        
        
    
    
    
