function DisplayPosteriorTracks(X,posteriorProb)
%-------------------------------------------------------------------------- 
% This function plots all of the particle tracks with color proportional to 
% posterior probability for each diffusive state
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

[numTracks numStates] = size(posteriorProb);

% set up a heat-based color map
bin = 100;
edges = 0:1/bin:1;
colorSet = hot(bin+1);
colorSet = colorSet(end:-1:1,:);

% figure('position',[50 50 1200 300]);  
for k = 1:numStates
    figure; hold on; box on;
    for p = 1:numTracks
        % find index of colormap which corresponds to posterior prob.
        [MIN index] = min((posteriorProb(p,k) - edges).^2);    
        plot(X{p}(:,1),X{p}(:,2),'color',colorSet(index,:));
    end
    axis tight;
    xlabel('Positions (\mum)','fontsize',16);
    ylabel('Positions (\mum)','fontsize',16);
    title(['State ' num2str(k)],'fontsize',16);
end
