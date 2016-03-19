%--------------------------------------------------------------------------
% This script calculates the expectation step: the posterior probabilities
% and the log-likelihood given the model parameters
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%--------------------------------------------------------------------------

clear all;
clc;
close all;
addpath('pEMv2');

%%  load file

[filename,dirpath] = uigetfile('*.mat','Select protein track positions mat file');
data = load(fullfile(dirpath,filename));
Xraw = data.X;

%% user set parameters

% movie parameters
dt = .032;              % time between steps
dE = .032;              % exposure time

% pEM parameters
minStates = 3;          % minimum number of states to explore
maxStates = 5;          % maximum number of states to explore
numReinitialize = 3;    % number of reinitialization trials
numPerturb = 5;         % number of perturbation trials
maxiter = 10000;        % maximum number of iterations within EM trial
convergence = 1e-7;     % convergence criteria for change in log-likelihood 
lambda = 0.00;          % shrinkage factor (useful when numerical issues calculating
                        % inverse of covariance matrix, labmda = 0.0 for no correction 
                        % lambda = 0.01 for correction)

numFeatures = 3;        % number of covariance features to include (min=2 for
                        % normal diffusion, 3-5 for non-normal diffusion)
splitLength = 20;       % length of steps to split each track

%% run pEM version 2

% split tracks into equal bin sizes
[X,splitIndex] = SplitTracks(Xraw,splitLength);

% structure for track info
trackInfo.numberOfTracks = length(X);   % number of tracks
trackInfo.dimensions = size(X{1},2);    % particle track dimensions
trackInfo.numFeatures = numFeatures;
trackInfo.splitLength = splitLength;
trackInfo.dt = dt;                      % frame duration
trackInfo.R = 1/6*dE/dt;                % motion blur coefficient
trackInfo.lambda = lambda;

% structure for pEM
params.minStates = minStates;
params.maxStates = maxStates;
params.numFeatures = numFeatures;
params.numReinitialize = numReinitialize;
params.numPerturbation = numPerturb;    % number of perturbations trials
params.converged = convergence;         % convergence condition for EM
params.maxiter = maxiter;               % maximum number of iterations for EM
params.verbose = 1;                     % display progress on command window (0,1)

% run pEMv2 
results = pEMv2_SPT(X,trackInfo,params); 

% display results
optimalSize = results.optimalSize;
optimalVacf = results.optimalVacf;
optimalP = results.optimalP;
disp('-------------------------------------------------------');
disp(['OptimalSize: ' num2str(optimalSize) ' states']);
for i = 1:numFeatures
    disp(['Sigma_k(i,i+' num2str(i-1) '): ' num2str(optimalVacf(:,i)') ' um^2']);
end
disp(['pi_k: ' num2str(optimalP)]);
disp('-------------------------------------------------------');

% save results
saveFolder = 'Results';
if ~isdir(saveFolder)
    mkdir(saveFolder)
end
[tmp, name] = fileparts(filename);
disp(['Saving results: Results/' name '.mat']); 
save(fullfile(saveFolder,[name '.mat']),'results');

%%











