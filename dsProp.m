

load neuronTable.mat
 
% Remove all the BAD cells
neuronTable = neuronTable(strcmp(neuronTable.cellID,"ON-OFF") | strcmp(neuronTable.cellID,"ON"),:);
neuronTable = neuronTable(strcmp(neuronTable.calciumSensor,"calDye"),:); %Only cal dye retinas, no gcamps

% Set a DSI or VARIANCE thresholds??
DSIthresh = 0; %Set to 0 for no Thresh
DSIsigThresh = 0.95; 
VARthresh = 10000; %Set to 10 for no Thresh; Set to 0.2 for reasonable thresh
dFoFthresh = 0;
neuronTable = neuronTable(neuronTable.DSI > DSIthresh & neuronTable.DSIsig > DSIsigThresh &...
    neuronTable.varSum < VARthresh & max(neuronTable.meanRespToBars,[],2) > dFoFthresh,:);

% Filter for just one retina location
neuronTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);

%In case you just want to look at one eye:
% neuronTable = neuronTable(strcmp(neuronTable.eye,"right"),:);



%% Loading the different groups

%ADULTS:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190116") | strcmp(neuronTable.animalID,"190117") |...
%     strcmp(neuronTable.animalID,"190404") | strcmp(neuronTable.animalID,"190405") |...
%     strcmp(neuronTable.animalID,"190419"),:); %The adults

adults = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"NR"),:);



%P14s:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190522") | strcmp(neuronTable.animalID,"190523"),:); %The devs
% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190523"),:); 

p14s = neuronTable(neuronTable.age < 17 & strcmp(neuronTable.condition,"NR"),:);




%DARK REARS:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190424") | strcmp(neuronTable.animalID,"190425") |...
%     strcmp(neuronTable.animalID,"190426") | strcmp(neuronTable.animalID,"190323")...
%     ,:); %The dark-rears

DRs = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"DR"),:);

%%




%% Identify the 4 DS clusters





% % Do you just want to look at HB9s?
% neuronTable = neuronTable(strcmp(neuronTable.GFPid, "drd4") ,:);
% neuronTable = neuronTable(strcmp(neuronTable.GFPid, "hb9") ,:);
% neuronTable = neuronTable(strcmp(neuronTable.cellID,"ON-OFF"),:);
% neuronTable = neuronTable(strcmp(neuronTable.cellID,"ON"),:);


[adults.x,adults.y,adults.idxDS, adults_C] = classDS(adults.prefDirCorr,adults.vecSum, 4);
[DRs.x,DRs.y,DRs.idxDS, DRs_C] = classDS(DRs.prefDirCorr,DRs.vecSum, 4);
[p14s.x,p14s.y,p14s.idxDS, p14s_C] = classDS(p14s.prefDirCorr,p14s.vecSum, 4);

scatterAxisLim = 0.8;

figure
subplot(1,3,1)
scatter(adults.x,adults.y,[],adults.idxDS)
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])
subplot(1,3,2)
scatter(adults.x(strcmp(adults.cellID,'ON-OFF')),adults.y(strcmp(adults.cellID,'ON-OFF')),[],adults.idxDS(strcmp(adults.cellID,'ON-OFF')))
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])
subplot(1,3,3)
scatter(adults.x(strcmp(adults.cellID,'ON')),adults.y(strcmp(adults.cellID,'ON')),[],adults.idxDS(strcmp(adults.cellID,'ON')))
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])

figure
subplot(1,3,1)
scatter(DRs.x,DRs.y,[],DRs.idxDS)
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])
subplot(1,3,2)
scatter(DRs.x(strcmp(DRs.cellID,'ON-OFF')),DRs.y(strcmp(DRs.cellID,'ON-OFF')),[],DRs.idxDS(strcmp(DRs.cellID,'ON-OFF')))
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])
subplot(1,3,3)
scatter(DRs.x(strcmp(DRs.cellID,'ON')),DRs.y(strcmp(DRs.cellID,'ON')),[],DRs.idxDS(strcmp(DRs.cellID,'ON')))
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])

figure
subplot(1,3,1)
scatter(p14s.x,p14s.y,[],p14s.idxDS)
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])
subplot(1,3,2)
scatter(p14s.x(strcmp(p14s.cellID,'ON-OFF')),p14s.y(strcmp(p14s.cellID,'ON-OFF')),[],p14s.idxDS(strcmp(p14s.cellID,'ON-OFF')))
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])
subplot(1,3,3)
scatter(p14s.x(strcmp(p14s.cellID,'ON')),p14s.y(strcmp(p14s.cellID,'ON')),[],p14s.idxDS(strcmp(p14s.cellID,'ON')))
axis([-scatterAxisLim scatterAxisLim -scatterAxisLim scatterAxisLim])


figure, scatter(p14s.x,p14s.y,[],p14s.idxDS)


propBars(adults)
propBars(DRs)
propBars(p14s)

%% find num cells per directions
% 1 2 3 4 =  T D N V

[adultCellNumPerDirPerFOV] = findNumCellsPerDirPerFOV(adults);
[p14CellNumPerDirPerFOV] = findNumCellsPerDirPerFOV(p14s);
[DRCellNumPerDirPerFOV] = findNumCellsPerDirPerFOV(DRs);

propAdultCellNumPerDirPerFOV = adultCellNumPerDirPerFOV ./ sum(adultCellNumPerDirPerFOV,2);
propP14CellNumPerDirPerFOV = p14CellNumPerDirPerFOV ./ sum(p14CellNumPerDirPerFOV,2);
propDRCellNumPerDirPerFOV = DRCellNumPerDirPerFOV ./ sum(DRCellNumPerDirPerFOV,2);

propAdultCellNumPerDirPerFOV = [propAdultCellNumPerDirPerFOV, ones(length(propAdultCellNumPerDirPerFOV(:,1)),1)];
propP14CellNumPerDirPerFOV = [propP14CellNumPerDirPerFOV, 2*ones(length(propP14CellNumPerDirPerFOV(:,1)),1)];
propDRCellNumPerDirPerFOV = [propDRCellNumPerDirPerFOV, 3*ones(length(propDRCellNumPerDirPerFOV(:,1)),1)];

figure, plotSpread({propAdultCellNumPerDirPerFOV(:,1),propP14CellNumPerDirPerFOV(:,1),propDRCellNumPerDirPerFOV(:,1)})

%% Functions

%Classify DS clusters based on prefDir

function [x,y,idxOrdered, C] = classDS(prefDir,vecSum, numClusters)

    [x,y] = pol2cart(prefDir,vecSum);
    [xNorm,yNorm] = pol2cart(prefDir,ones(length(prefDir),1));

    [idx,C] = kmeans([xNorm yNorm],numClusters);

    idxOrdered = nan(length(idx),1);

    % Reorder grps based on max X or Y
    %Temporal = 1
    [dummy,id] = max(C(:,1));
    idxOrdered(idx == id) = 1;
    %Dorsal = 2
    [dummy,id] = max(C(:,2));
    idxOrdered(idx == id) = 2;
    %Nasal = 3
    [dummy,id] = min(C(:,1));
    idxOrdered(idx == id) = 3;
    %Ventral = 4
    [dummy,id] = min(C(:,2));
    idxOrdered(idx == id) = 4;
    %Done with Reorder

end



%Build a proportion bar plot of number of DS cells in each direction

function propBars(tempTable)

    props = [length(find(tempTable.idxDS == 1)) length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON-OFF'))) length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON')));...
        length(find(tempTable.idxDS == 2)) length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON-OFF'))) length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON')));...
        length(find(tempTable.idxDS == 3)) length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON-OFF'))) length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON')));...
        length(find(tempTable.idxDS == 4)) length(find(tempTable.idxDS == 4 & strcmp(tempTable.cellID,'ON-OFF'))) length(find(tempTable.idxDS == 4 & strcmp(tempTable.cellID,'ON')))]...
        ./[length(tempTable.idxDS) length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF'))) length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))];
    
    figure
    bar(props','stacked', 'FaceColor', 'none','BarWidth', 0.5)
    
    % Adding text next to bars that tell you number of DS cells
    text(0.5,length(find(tempTable.idxDS == 1))/length(tempTable.idxDS)/2,num2str(length(find(tempTable.idxDS == 1))))
    text(0.5,length(find(tempTable.idxDS == 1))/length(tempTable.idxDS)+length(find(tempTable.idxDS == 2))/length(tempTable.idxDS)/2,num2str(length(find(tempTable.idxDS == 2))))
    text(0.5,length(find(tempTable.idxDS == 1))/length(tempTable.idxDS)+length(find(tempTable.idxDS == 2))/length(tempTable.idxDS)+length(find(tempTable.idxDS == 3))/length(tempTable.idxDS)/2,num2str(length(find(tempTable.idxDS == 3))))
    text(0.5,length(find(tempTable.idxDS == 1))/length(tempTable.idxDS)+length(find(tempTable.idxDS == 2))/length(tempTable.idxDS)+length(find(tempTable.idxDS == 3))/length(tempTable.idxDS)+length(find(tempTable.idxDS == 4))/length(tempTable.idxDS)/2,num2str(length(find(tempTable.idxDS == 4))))
    text(1.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))/2,num2str(length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON-OFF')))))
    text(1.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))+length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))/2,num2str(length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON-OFF')))))
    text(1.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))+length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))+length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))/2,num2str(length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON-OFF')))))
    text(1.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))+length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))+length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))+length(find(tempTable.idxDS == 4 & strcmp(tempTable.cellID,'ON-OFF')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))/2,num2str(length(find(tempTable.idxDS == 4 & strcmp(tempTable.cellID,'ON-OFF')))))
    text(2.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))/2,num2str(length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON')))))
    text(2.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))+length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))/2,num2str(length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON')))))
    text(2.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))+length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))+length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))/2,num2str(length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON')))))
    text(2.5,length(find(tempTable.idxDS == 1 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))+length(find(tempTable.idxDS == 2 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))+length(find(tempTable.idxDS == 3 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))+length(find(tempTable.idxDS == 4 & strcmp(tempTable.cellID,'ON')))/length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))/2,num2str(length(find(tempTable.idxDS == 4 & strcmp(tempTable.cellID,'ON')))))
    
    text(0.5,0.01,num2str(length(tempTable.idxDS)))
    text(1.5,0.01,num2str(length(tempTable.idxDS(strcmp(tempTable.cellID,'ON-OFF')))))
    text(2.5,0.01,num2str(length(tempTable.idxDS(strcmp(tempTable.cellID,'ON')))))
    
end


% proportion of temporal cells per animal

function [cellNumbPerDirPerFOV] = findNumCellsPerDirPerFOV(tempTable)

listIDs = unique(tempTable.animalID);
cellNumbPerDirPerFOV = zeros(length(listIDs),4); % T D N V
for i = 1:length(listIDs)
    funTempTable = tempTable(strcmp(tempTable.animalID, listIDs(i)),:); 
    
    cellNumbPerDirPerFOV(i,1) = length(funTempTable.idxDS(funTempTable.idxDS == 1));
    cellNumbPerDirPerFOV(i,2) = length(funTempTable.idxDS(funTempTable.idxDS == 2));
    cellNumbPerDirPerFOV(i,3) = length(funTempTable.idxDS(funTempTable.idxDS == 3));
    cellNumbPerDirPerFOV(i,4) = length(funTempTable.idxDS(funTempTable.idxDS == 4));
    
end

end





