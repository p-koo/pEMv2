clear all;
clc;
close all;

%%

[filename,dirpath] = uigetfile('*.mat');
load(fullfile(dirpath,filename));

%%

offset = .05;
numStates = 4;
colors = 'bgrc';
binIndex = results.binIndex;
gamma = results.state(numStates).pem.gamma;
binSize = size(results.binX{1},1);

for i = 1:length(results.X)
    X = results.X{i};

    % classified states
    index = find(binIndex == i);
    state = mean(gamma(index,:,:),3);
    [MAX,classifySeq] = max(state,[],2);     
    classSeq = [classifySeq(1)];
    for j = 1:length(index)
        range = j*binSize-binSize+1:j*binSize;
        classSeq = [classSeq; ones(binSize,1)*classifySeq(j)];
    end
    N = length(classSeq);
    
    % true states
    trueSeq = results.stateSeq{i};
    trueSeq = trueSeq(1:N);
    
    X = X(1:N,:);
    minx = min(X);
    maxx = max(X);
    
    figure(1); clf; 
    subplot(1,2,1); box on; hold on; axis equal;
    for j = 1:N-1
        plot(X(j:j+1,1),X(j:j+1,2),'color',colors(trueSeq(j)));
    end
    axis([minx(1)-offset maxx(1)+offset minx(2)-offset maxx(2)+offset]);
    set(gca,'linewidth',1.5,'xticklabel',[],'yticklabel',[]);
    
    subplot(1,2,2); box on; hold on; axis equal;
    for j = 1:N-1
        plot(X(j:j+1,1),X(j:j+1,2),'color',colors(classSeq(j)));
    end    
    axis([minx(1)-offset maxx(1)+offset minx(2)-offset maxx(2)+offset]);
    set(gca,'linewidth',1.5,'xticklabel',[],'yticklabel',[]);
    rectangle('position',[minx(1) minx(2) .2 .02],'facecolor','k','edgecolor','k');
    title(num2str(i));
    
    pause;
    
end


%%

plotN = 112;

dirpath = '';
filename{1} = 'Results_Results_VaryTransition_case3_7_1_5';
filename{2} = 'Results_Results_VaryTransition_case3_7_1_15';
filename{3} = 'Results_Results_VaryTransition_case3_7_1_30';

trackIndex = [14 61];

for z = 1:length(filename)
    load(fullfile(dirpath,filename{z}));

    offset = .05;
    numStates = 4;
    colors = 'bgrc';
    binIndex = results.binIndex;
    gamma = results.state(numStates).pem.gamma;
    binSize = size(results.binX{1},1)+1;

    for i = [14 61]
        X = results.X{i};

        % classified states
        index = find(binIndex == i);
        state = mean(gamma(index,:,:),3);
        [MAX,classifySeq] = max(state,[],2);     
         if z == 3
            index1 = find(classifySeq == 1);
            index2 = find(classifySeq == 2);
            classifySeq(index1) = 2;
            classifySeq(index2) = 1;
        end
        classSeq = [classifySeq(1)];
        for j = 1:length(index)
            range = j*binSize-binSize+1:j*binSize;
            classSeq = [classSeq; ones(binSize,1)*classifySeq(j)];
        end
        classSeq(end) = [];
        N = length(classSeq);

        % true states
        trueSeq = results.stateSeq{i};
        trueSeq = trueSeq(1:N);

        X = X(1:plotN,:);
        minx = min(X);
        maxx = max(X);

        figure(1); clf; 
        subplot(1,2,1); box on; hold on; axis equal;
        for j = 1:plotN-1
            plot(X(j:j+1,1),X(j:j+1,2),'color',colors(trueSeq(j)));
        end
        axis([minx(1)-offset maxx(1)+offset minx(2)-offset maxx(2)+offset]);
        rectangle('position',[minx(1) minx(2) .2 .02],'facecolor','k','edgecolor','k');
        set(gca,'linewidth',1.5,'xticklabel',[],'yticklabel',[]);

        subplot(1,2,2); box on; hold on; axis equal;
        for j = 1:plotN-1
            plot(X(j:j+1,1),X(j:j+1,2),'color',colors(classSeq(j)));
        end    
        axis([minx(1)-offset maxx(1)+offset minx(2)-offset maxx(2)+offset]);
        set(gca,'linewidth',1.5,'xticklabel',[],'yticklabel',[]);
        print('-depsc',[filename{z} '_X_' num2str(i)]);

        figure(2); clf; hold on; box on;
        plot(1:plotN,trueSeq(1:plotN),'b','linewidth',3);
        plot(1:plotN,classSeq(1:plotN),':r','linewidth',3);
        set(gca,'linewidth',2,'fontsize',35);
        axis([1 plotN .9 4.1]);
        set(gca,'ytick',[1 2 3 4]);
        print('-depsc',[filename{z} '_state_' num2str(i)]);

    end
end


%%

dirpath = '';
filename{1} = 'Results_Results_VaryTransition_case3_7_1_5';
filename{2} = 'Results_Results_VaryTransition_case3_7_1_15';
filename{3} = 'Results_Results_VaryTransition_case3_7_1_30';

trackIndex = [1 14 20 36 61 75 76 88 101 109];

offset = .05;
numStates = 4;
colors = 'grc';
style = {':',':',':'};

for i = [1 14 20 36 61 75 76 88 101 109]
    figure(2); clf; hold on; box on;
    z = 3;
    load(fullfile(dirpath,filename{z}));
    X = results.X{i};
    binIndex = results.binIndex;
    gamma = results.state(numStates).pem.gamma;
    binSize = size(results.binX{1},1);

    % classified states
    index = find(binIndex == i);
    state = mean(gamma(index,:,:),3);
    [MAX,classifySeq] = max(state,[],2);     
    if z == 3
        index1 = find(classifySeq == 1);
        index2 = find(classifySeq == 2);
        classifySeq(index1) = 2;
        classifySeq(index2) = 1;
    end

    classSeq = [classifySeq(1)];
    for j = 1:length(index)
        range = j*binSize-binSize+1:j*binSize;
        classSeq = [classSeq; ones(binSize,1)*classifySeq(j)];
    end
    N = length(classSeq);

    % true states
    trueSeq = results.stateSeq{i};
    trueSeq = trueSeq(1:N);

    plot(1:N,trueSeq,'b','linewidth',2);

    for z = 1:length(filename)
        load(fullfile(dirpath,filename{z}));
        X = results.X{i};
        binIndex = results.binIndex;
        gamma = results.state(numStates).pem.gamma;
        binSize = size(results.binX{1},1);

        % classified states
        index = find(binIndex == i);
        state = mean(gamma(index,:,:),3);
        [MAX,classifySeq] = max(state,[],2);     
        if z == 3
            index1 = find(classifySeq == 1);
            index2 = find(classifySeq == 2);
            classifySeq(index1) = 2;
            classifySeq(index2) = 1;
        end
        
        classSeq = [classifySeq(1)];
        for j = 1:length(index)
            range = j*binSize-binSize+1:j*binSize;
            classSeq = [classSeq; ones(binSize,1)*classifySeq(j)];
        end
        N = length(classSeq);

        % true states
        trueSeq = results.stateSeq{i};
        trueSeq = trueSeq(1:N);

        X = X(1:N,:);
        minx = min(X);
        maxx = max(X);

%         if z == 1
%         end
        plot(1:N,classSeq,[style{z} colors(z)],'linewidth',2);
    end

    
    set(gca,'linewidth',2,'fontsize',25);
    axis([1 N .9 4.1]);
    set(gca,'ytick',[1 2 3 4]);
    print('-depsc',['State_' num2str(i)]);
end






