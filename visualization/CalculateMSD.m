function msd = CalculateMSD(x,numLags)

dim = size(x,2);

msd = zeros(1,numLags);
for i = 1:numLags
    rho = [];
    for j = 1:i
        for k = 1:dim
            rho = [rho; diff(x(j:i:end,k)).^2];
        end
    end
    msd(i) = mean(rho);
end

