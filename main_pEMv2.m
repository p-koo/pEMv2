%--------------------------------------------------------------------------
% This script runs perturbation expectation-maximization version 2 (pEMv2) 
% on a set of simulated particle tracks. The tracks have to be stored in a 
% mat file under the variable X, which is a cell that contains each
% trajectory X{1} = [x_1 y_1], X{2} = [x_2, y_2]... where x_i and y_i are
% vectors of the positions of the trajectories.  The output of pEMv2 is
% saved in a mat file in the results folder.
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
X_raw = data.X;

% get only one dimensional tracks
X_new = cell(length(X_raw), 1);
for i=1:length(X_raw)
    X_new{i} = X_raw{i}(:,1);
end
X = X_new;


%% user set parameters

% movie parameters
dt = .032;              % time between steps
dE = .032;              % exposure time

% pEM parameters
minStates = 2;          % minimum number of states to explore
maxStates = 6;          % maximum number of states to explore
numReinitialize = 3;    % number of reinitialization trials
numPerturb = 30;        % number of perturbation trials
maxiter = 10000;        % maximum number of iterations within EM trial
convergence = 1e-5;     % convergence criteria for change in log-likelihood 
lambda = 0.00;          % shrinkage factor (useful when numerical issues calculating
                        % inverse of covariance matrix, labmda = 0.0 for no correction 
                        % lambda = 0.01 for correction)


splitLength = 15;                   % length of steps to split each track
numFeatures = 6; %splitLength-1;        % number of covariance features to include (min=2 for
                                    % normal diffusion, 3-5 for non-normal diffusion)
%% run pEM version 2

% split tracks into equal bin sizes
[X,splitIndex] = SplitTracks(X_raw,splitLength);

% structure for track info
trackInfo.numberOfTracks = length(X);   % number of tracks
trackInfo.dimensions = size(X{1},2);    % particle track dimensions
trackInfo.numFeatures = numFeatures;    % number of features to retain in covariance matrix
trackInfo.splitLength = splitLength;    % length of each bin
trackInfo.splitIndex = splitIndex;      % index of each track
trackInfo.dt = dt;                      % frame duration
trackInfo.R = 1/6*dE/dt;                % motion blur coefficient
trackInfo.lambda = lambda;              % shrinkage factor


% structure for pEM
params.minStates = minStates;               % minimum number of states to try
params.maxStates = maxStates;               % maximum number of states to try
params.numFeatures = numFeatures;           % number of features in covariance elements
params.numReinitialize = numReinitialize;   % number of reinitialization trials
params.numPerturbation = numPerturb;        % number of perturbations trials
params.converged = convergence;             % convergence condition for EM
params.maxiter = maxiter;                   % maximum number of iterations for EM
params.verbose = 1;                         % display progress on command window (0,1)


% calculate the displacements for each particle track
deltaX = cell(trackInfo.numberOfTracks,1);
for i = 1:trackInfo.numberOfTracks
    deltaX{i} = diff(X{i});
end

% calculate relevant properties to enhance compuatational time
[trackInfo.vacf_exp,trackInfo.xbar_exp] = CovarianceProperties(deltaX,numFeatures);

% run pEMv2 
results = pEMv2_SPT(deltaX,trackInfo,params); 

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

% save results in Results/filename
[tmp, name] = fileparts(filename);
saveFolder = fullfile('results',name);
if ~isdir(saveFolder)
    mkdir(saveFolder)
end
disp(['Saving results: ' fullfile(saveFolder,['results.mat'])]); 
save(fullfile(saveFolder,['results.mat']),'results');


%%  Fit covariance matrices with theoretical models 

addpath('CovarianceModels');

gamma= results.posteriorProb;
[argvalue, argmax] = max(gamma');

params.dt = dt;
params.R = 1/6*dE/dt;

method = 'other';
models = {'normal', 'confined', 'fBM'};
states = struct;
for j = 1:optimalSize
    disp(['state ' int2str(j)])
    C_exp = toeplitz(optimalVacf(j,:));
    params.D0 = (C_exp(1,1) + 2*C_exp(1,2))/(2*dt);
    params.sigma0 = sqrt(abs(C_exp(1,1)/2 - params.D0*dt*(1-2*params.R)));
    params.L0 = params.D0*5;
    params.A0 = 1;
    states(j).fit = struct;
    for i = 1:length(models)
        [p, loss] = OptimalParameters(C_exp, models{i}, params);
        states(j).fit(i).name = models(i);
        states(j).fit(i).p = p;
        states(j).fit(i).loss = loss;
        states(j).fit(i).N = length(p);
        disp(['    model: ' models{i}]);
        disp(['         p: ' num2str(p)]);
        disp(['         loss: ' num2str(loss)]);
    end
end


%%
