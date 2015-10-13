function [C2] = shrinkcov(C,lambda)
% Ledoit-Wolf optimal shrinkage estimator for cov(X)  

% shrunk final estimate C
C2 = (1-lambda)*C + lambda*diag(diag(C));
        
end
