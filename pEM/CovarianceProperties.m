function [vacf_exp,xbar_exp] = CovarianceProperties(deltaX,numFeatures)
%-------------------------------------------------------------------------- 
% Summary: This function calculates the empirical covariance properties for
% each particle track, assuming normal diffusion with a constant
% diffusivity and static localization noise throughout each particle track.
% Each dimension is calculated separately.
% 
% Input:
%       deltaX = cell with the particle track displacements
%
% Output:
%       diagonals = average diagonal elements of empirical covariance matrix
%       correlations = average nearest-neighbor correlations of empirical
%       covariance matrix for each particle track and each dimension
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

numTracks = length(deltaX);
dim = size(deltaX{1},2);

% calculate velocity auto-covariance for all particle tracks
vacf_exp = zeros(numTracks,numFeatures,dim); 
for j = 1:dim
    for i = 1:numTracks
        diffX = deltaX{i}(:,j);            
        C = diffX*diffX';
        tmpvacf = CalculateVACF(C);
        vacf_exp(i,:,j) = tmpvacf(1:numFeatures);
    end
end 

xbar_exp = zeros(numTracks,dim);
for i = 1:numTracks
    xbar_exp(i,:) = mean(deltaX{i});
end


