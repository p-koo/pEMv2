function [meanmsd,errormsd,nmsd] = AverageMSD(truthX,numStates,numLags)

meanmsd = zeros(numStates,numLags);
for p = 1:numStates
    X = truthX(p).X;    
    msd = zeros(length(X),numLags);
    for k = 1:length(X)
        x = X{k}(:,1);
        y = X{k}(:,2);

        if length(x) > 1
            for i = 1:numLags
                rho = [];
                for j = 1:i
                    rho = [rho; diff(x(j:i:end)).^2 + diff(y(j:i:end)).^2];
                end
                msd(k,i) = mean(rho);
            end
        end
    end
    
    for k = 1:numLags
        index = find(~isnan(msd(:,k)));
        meanmsd(p,k) = mean(msd(index,k));
        errormsd(p,k) = std(msd(index,k));
        nmsd(p,k) = length(index);
    end
end



