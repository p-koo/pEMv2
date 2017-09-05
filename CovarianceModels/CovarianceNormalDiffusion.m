function C = CovarianceNormalDiffusion(N,D,sigma,R,dt)

% generate covariance matrix
msd = 2*D*dt;
staticnoise = 2*sigma^2;
motionblur = -2*R*msd;
C = toeplitz([msd+staticnoise+motionblur -(staticnoise+motionblur)/2 zeros(1,N-2)]);
