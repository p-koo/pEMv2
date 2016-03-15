clear all;
clc; 
close all;

%%

[filename,dirpath] = uigetfile('*.mat','Select file');
load(fullfile(dirpath,filename));

numStates = 4;

X = results.X;
numTracks = length(X);

% true state
state = zeros(numTracks,1);
for i = 1:numTracks
    state(i) = results.stateSeq{i}(1);
end

% classified state
posteriorProb = mean(results.state(numStates).posteriorProb,3);
[MAX,classify] = max(posteriorProb,[],2);

trueState = zeros(numTracks,1);
classifyState = zeros(numTracks,1);
k = 1;
for i = 1:numStates
    index = find(state == i);
    trueState(k:k+length(index)-1) = state(index);
    classifyState(k:k+length(index)-1) = classify(index);
    k = k + length(index);
end
sum(trueState==classifyState)/numTracks

%%
space = .6;
width = 30;
height = 50;

colorSet = 'rgbcm';

figure; hold on; box on;
k = 1;
for i = 1:height
    for j = 1:width
        x = X{k}(:,1) + i*space;
        y = X{k}(:,2) + j*space;
        plot(x,y,colorSet(trueState(k)));
        k = k + 1;
    end
end
axis equal;
axis([-1 height*space+1 -1.5 width*space+1]);
set(gca,'xtick',[],'ytick',[]);
rectangle('position',[-.5 -1 5 .3],'facecolor','k');
print('-depsc','case1_tracks_true');

figure; hold on; box on;
k = 1;
for i = 1:height
    for j = 1:width
        x = X{k}(:,1) + i*space;
        y = X{k}(:,2) + j*space;
        plot(x,y,colorSet(classifyState(k)));
        k = k + 1;
    end
end
axis equal;
axis([-1 height*space+1 -1.5 width*space+1]);
set(gca,'xtick',[],'ytick',[]);
print('-depsc','case1_tracks_classify');


%%









        
        