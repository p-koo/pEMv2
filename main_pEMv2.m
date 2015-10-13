clear all;
clc;
close all;
addpath('pEM');

%%  load file

[filename dirpath] = uigetfile('*.mat','Select protein track positions mat file');
data = load(fullfile(dirpath,filename));
Xraw = data.X;


%% user set parameters

% movie parameters
dt = .032;
dE = .032;

% pEM parameters
minStates = 3;          % minimum number of states to explore
maxStates = 5;          % maximum number of states to explore
numReinitialize = 3;
numPerturb = 5;        % number of perturbation trials
maxiter = 10000;        % maximum number of iterations within EM trial
convergence = 1e-7;     % convergence criteria for change in log-likelihood 
lambda = 0.00;

numFeatures = 2;
splitLength = 10;

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

results = pEM_SPT(X,trackInfo,params); 

%% display results

optimalSize = results.optimalSize;
optimalVacf = results.optimalVacf;
optimalP = results.optimalP;
disp('-------------------------------------------------------');
disp(['OptimalSize: ' num2str(optimalSize)]);
disp(['D_k: ' num2str((optimalVacf(:,1) + 2*optimalVacf(:,2))'/2/dt)]);
disp(['pi_k: ' num2str(optimalP)]);
disp('-------------------------------------------------------');


%%












