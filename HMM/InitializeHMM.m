function hmmmodel = InitializeHMM(deltaX, pik, basevacf, prob)

[K,T] = size(basevacf);
D = length(deltaX{1});

% initial state
p = repmat(pik',1,length(deltaX));

% transition matrix
a = InitializeTransitionMatrix(K,prob);
    
% initial sigma
sigma = zeros(D,D,K);
for i = 1:K
    sigma(:,:,i) = toeplitz([basevacf(i,:) zeros(1,D-T)]);
end

% Emission probabilities
b = EmissionProbabilities(deltaX, sigma, 0);

hmmmodel.p = p;
hmmmodel.a = a;
hmmmodel.b = b;
hmmmodel.sigma = sigma;

