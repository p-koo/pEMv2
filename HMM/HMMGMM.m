function hmmmodel = HMMGMM(splitX, hmmmodel, trackInfo, params)

maxiter = params.maxiter;
tolerance = params.tolerance;
verbose = params.verbose;
splitIndex = trackInfo.splitIndex;

p = hmmmodel.p;
a = hmmmodel.a;
b = hmmmodel.b;
sigma = hmmmodel.sigma;

K = size(a,1);

pbest = p;
abest = a;
bbest = b;
sigmabest = sigma;
gammankbest = 0;

% HMM model parameters
Lold = -Inf;
for i = 1:maxiter
    
    % expectation step
    [est,Lnew,gammank] = ExpectationHMM(p,a,b,splitIndex);
    disp([num2str(i) ': ' num2str(Lnew)]);

    % maximization step
    [p,a,b,sigma] = MaximizationHMM(splitX,gammank,est,K,trackInfo);
    
    % check percentage similar
    [MAX,index] = max(gammank,[],2);

    if (Lnew - Lold)/abs(Lold) < tolerance
        if verbose == 1
            disp('Converged...');
            disp(['iteration ' num2str(i) ': logL = ' num2str(Lnew)]);
        end
        break;
    else
        Lold = Lnew;
        pbest = p;
        abest = a;
        bbest = b;
        sigmabest = sigma;
        gammankbest = gammank;
    end
end

if i == maxiter
    if verbose == 1
        disp('Did not converge...');
    end
end

hmmmodel.p = pbest;
hmmmodel.a = abest;
hmmmodel.b = bbest;
hmmmodel.sigma = sigmabest;
hmmmodel.gammank = gammankbest;
hmmmodel.logL = Lold;

