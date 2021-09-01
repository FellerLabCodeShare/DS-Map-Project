
load ventralAllClusteredWithLocs.mat
load b2koClass.mat

b2koClass = b2koClass(b2koClass.age > 29,:);

tempTable = NR;



% Set a DSI or VARIANCE thresholds??
DSIthresh = 0; %Set to 0 for no Thresh
DSIsigThresh = 0.95; 
VARthresh = 10000; %Set to 10 for no Thresh; Set to 0.2 for reasonable thresh
dFoFthresh = 0;
tempTable = tempTable(tempTable.DSI > DSIthresh & tempTable.DSIsig > DSIsigThresh &...
    tempTable.varSum < VARthresh & max(tempTable.meanRespToBars,[],2) > dFoFthresh,:);

tempTable.angleID = deg2rad(tempTable.idxDS*90 - 90);
tempTable.deltaTheta = wrapToPi(tempTable.prefDirCorr - tempTable.angleID);
% tempTable.deltaTheta = wrapToPi(tempTable.prefDir - tempTable.angleID);

listFOVs = unique(tempTable.fileName);

meanDeviations = nan(length(listFOVs),4);
distFromON = nan(length(listFOVs),1);


for i = 1:length(listFOVs)
   fovTable = tempTable(strcmp(tempTable.fileName,listFOVs(i)),:);
   meanDeviations(i,1) = mean(fovTable.deltaTheta(fovTable.idxDS == 1));
   meanDeviations(i,2) = mean(fovTable.deltaTheta(fovTable.idxDS == 2));
   meanDeviations(i,3) = mean(fovTable.deltaTheta(fovTable.idxDS == 3));
   meanDeviations(i,4) = mean(fovTable.deltaTheta(fovTable.idxDS == 4));
    
   if strcmp(fovTable.location,"ventroTemporal") == 1
       distFromON(i) = mean(fovTable.distanceFromON)*-1;
   else
       distFromON(i) = mean(fovTable.distanceFromON);
   end
   
end



figure('Name','angle deviations','NumberTitle','off')
set(gcf, 'Position', [10   10   900   900]);

subplot(2,2,1)
    hold
    scatterWithLine(distFromON,meanDeviations(:,1))
subplot(2,2,2)
    hold
    scatterWithLine(distFromON,meanDeviations(:,2))
subplot(2,2,3)
    hold
    scatterWithLine(distFromON,meanDeviations(:,3))
subplot(2,2,4)
    hold
    scatterWithLine(distFromON,meanDeviations(:,4))

%% Calc deviation in central and peripheral retina

central = -1000;
peripheral = 1000;

grpIds = nan(size(distFromON));
grpIds(distFromON<central) = 1;
grpIds(distFromON>peripheral) = 2;


figure('Name','angle deviations','NumberTitle','off')
set(gcf, 'Position', [10   10   900   450]);
subplot(1,4,1)
hold
boxplot(rad2deg(meanDeviations(:,1)),grpIds,'Notch','on')
ylim([-90 90]);
subplot(1,4,2)
hold
boxplot(rad2deg(meanDeviations(:,2)),grpIds,'Notch','on')
ylim([-90 90]);
subplot(1,4,3)
hold
boxplot(rad2deg(meanDeviations(:,3)),grpIds,'Notch','on')
ylim([-90 90]);
subplot(1,4,4)
hold
boxplot(rad2deg(meanDeviations(:,4)),grpIds,'Notch','on')
ylim([-90 90]);

% centralMeanDevs = meanDeviations(abs(distFromON) < central,:);
% peripheralMeanDevs = meanDeviations(abs(distFromON) > peripheral,:);


% figure('Name','angle deviations','NumberTitle','off')
% set(gcf, 'Position', [10   10   900   450]);
% subplot(1,2,1)
% hold
% plotSpread(abs(centralMeanDevs));
% subplot(1,2,2)
% hold
% plotSpread(abs(peripheralMeanDevs));

% rad2deg([nanmean(abs(centralMeanDevs)), nanmean(abs(peripheralMeanDevs))])



%% LinearCorrelation followed by linear regression to figure out if diff b/w exp grps in slopes



linReg = fitlm(distFromON,meanDeviations(:,2));



%% Functions



   function scatterWithLine(x,y)
    
    xWin = 2000;
   
    scatter(x,y,'k.')
    p = polyfit(x(~isnan(y)),y(~isnan(y)),1);
    x1 = linspace(-xWin,xWin);
    y1 = polyval(p,x1);
    plot(x1,y1,'k')
    
    ylim([-pi/2 pi/2]);
    xlim([-xWin xWin]);
    
    text(-xWin/2,pi/3,[num2str(1000*rad2deg(p(1))), ' deg/mm'])
   
   
   end
   
   
   
