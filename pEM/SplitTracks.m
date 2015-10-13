function [splitX,splitIndex] = SplitTracks(X,splitLength)

n = 1;
splitX = {};
splitIndex = [];
for i = 1:length(X)  
    N = length(X{i});
    k = 1;
    status = 1;
    while status == 1
        range = k*splitLength-splitLength+2-k:k*splitLength+1-k;
        
        if range(end) <= N    
            splitX{n} = X{i}(range,:);
            splitIndex = [splitIndex; i];
            k = k + 1;
            n = n + 1;
        else
            status = 0;
        end
    end        
end
