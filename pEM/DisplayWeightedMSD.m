function DisplayWeightedMSD(X,posteriorProb,numLags,dt)
%-------------------------------------------------------------------------- 
% This function plots the time averaged mean square displacement for all of 
% the particle tracks weighed by the posterior probability to be in a given 
% diffusive state
% 
% Code written by: 
%       Peter Koo
%       Yale University, Department of Physis, New Haven, CT, 06511  
%-------------------------------------------------------------------------- 

[numTracks numStates] = size(posteriorProb);

% calculate the time-averaged MSD
rho = zeros(numTracks,numLags);
for m = 1:numTracks
    for i = 1:numLags
        msd = [];
        for j = 1:i
            msd =  [msd; diff(X{m}(j:i:end,1)).^2 + diff(X{m}(j:i:end,2)).^2];
        end
        rho(m,i) = mean(msd);
    end
end

% plot posterior weighted taMSD
figure; hold on; box off;
colorSet = hsv(numStates);
timeLags = (1:numLags)*dt;
legendname = 'h = legend(';
for k = 1:numStates
    weight = posteriorProb(:,k)/sum(posteriorProb(:,k));
    weightedMSD = sum(rho.*repmat(weight,1,numLags));
    weightedMSDsquare = sum(rho.^2.*repmat(weight,1,numLags));
    errorMSD = sqrt(weightedMSDsquare - weightedMSD.^2);
    errorbar(timeLags,weightedMSD,errorMSD,'color',colorSet(k,:),'linewidth',1.1);
    legendname = [legendname '''State ' num2str(k) ''','];
end
legendname = [legendname(1:end-1) ');'];
eval(legendname);
set(gca,'fontsize',16,'linewidth',1.5);
set(h,'location','northwest','fontsize',14,'box','off');
xlabel('Time Lags (s)','fontsize',16);
ylabel('MSD','fontsize',16);

