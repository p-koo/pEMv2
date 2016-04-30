function b = EmissionProbabilities(deltaX,sigma,lambda)

T = length(deltaX);
[D,dim] = size(deltaX{1});
K = size(sigma,3);

% just use likelihood for normal diffusion
b = zeros(K,T);
for j = 1:K
    C = ShrinkCov(sigma(:,:,j), lambda);
    [logdetC,invC]= LogDeterminant(C);    
    for i = 1:dim
        for t = 1:T
            delta = deltaX{t}(:,i);
            L = -D/2*log(2*pi) -.5*logdetC - .5*delta'*invC*delta;
            b(j,t) = b(j,t) + L/dim; 
        end
    end
end

b = b - repmat(max(b),K,1);
b = exp(b)+eps;
b(isnan(b)) = 0;
b = b./(ones(K,1)*sum(b));




