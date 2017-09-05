function  [p, logL] = OptimalParameters(C_exp, model, params)
warning off;

% optimization parameters 'algorithm','sqp',
fitOptions = optimset('Display','off','TolX',1e-6,'MaxFunEvals',1000,'MaxIter',1000);

% particle track displacements
N = length(C_exp);

% initial parameters
dt = params.dt;
R = params.R;
D0 = params.D0;
sigma0 = params.sigma0;
L0 = params.L0;
A0 = params.A0;

% bounds for parameters
%Dmin = 1e-6;
%Dmax = 5;
%Smin = 1e-3;
%Smax = .5;
%Lmin = 1e-2;
%Lmax = 5;
%Amin = .1;
%Amax = 1.9;


switch model
    case 'normal'
        %lb = [Dmin Smin];
        %ub = [Dmax Smax];
        initialP0 = [D0 sigma0];
        %[p,logL,e,o,l,g,hessian] = fmincon(@(b) sum(sum((C_exp-CovarianceNormalDiffusion(N,b(1),b(2),R,dt)).^2)),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        [p,logL] = fminunc(@(b)sum(sum((C_exp-CovarianceNormalDiffusion(N,b(1),b(2),R,dt)).^2)),initialP0,fitOptions);
        p(2) = abs(p(2));
    case 'confined'
        %lb = [Dmin Lmin Smin];
        %ub = [Dmax Lmax Smax];
        initialP0 = [D0 L0 sigma0];
        %[p,logL,e,o,l,g,hessian] = fmincon(@(b)sum(sum((C_exp-CovarianceConfinedDiffusion(N,b(1),b(2),b(3),R,dt)).^2)),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        [p,logL] = fminunc(@(b)sum(sum((C_exp-CovarianceConfinedDiffusion(N,b(1),b(2),b(3),R,dt)).^2)),initialP0,fitOptions);
        p(3) = abs(p(3));
    case 'fBM'
        %lb = [Dmin Amin Smin];
        %ub = [Dmax Amax Smax];
        initialP0 = [D0 A0 sigma0];
        %[p,logL,e,o,l,g,hessian] = fmincon(@(b)sum(sum((C_exp-CovarianceFBM(N,b(1),b(2),b(3),R,dt)).^2)),initialP0,[],[],[],[],lb,ub,[],fitOptions);
        [p,logL] = fminunc(@(b)sum(sum((C_exp-CovarianceFBM(N,b(1),b(2),b(3),R,dt)).^2)),initialP0,fitOptions);
        p(3) = abs(p(3));
    otherwise
        disp('no case found');
end

