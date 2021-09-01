


load neuronTable.mat
load b2koClass.mat

b2koClass(strcmp(b2koClass.animalID,"201222"),:)=[];

% neuronTable = b2koClass; %If you want to analyze the B2s

% neuronTable = neuronTable(strcmp(neuronTable.GFPretina,"hb9"),:);
neuronTable = neuronTable(strcmp(neuronTable.cellID,"ON-OFF") | strcmp(neuronTable.cellID,"ON"),:);
neuronTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") | strcmp(neuronTable.location,"ventroNasal"),:);
neuronTable = neuronTable(strcmp(neuronTable.calciumSensor,"calDye")  | strcmp(neuronTable.calciumSensor,"calBryte"),:);

% neuronTable.location(:,1) = 'ventroNasal';
% 
% neuronTable = neuronTable(strcmp(neuronTable.animalID,"201101"),:);
% 
% neuronTable = neuronTable(strcmp(neuronTable.fileName,"190405_026_dFoF.mat") | strcmp(neuronTable.fileName,"190405_026_dFoF.mat") |...
%     strcmp(neuronTable.fileName,"190405_028_dFoF.mat") | strcmp(neuronTable.fileName,"190405_030_dFoF.mat") |...
%     strcmp(neuronTable.fileName,"190405_032_dFoF.mat") | strcmp(neuronTable.fileName,"190405_034_dFoF.mat"),:);

%In case you just want to look at one eye:
% neuronTable = neuronTable(strcmp(neuronTable.eye,"left"),:);

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190426") | strcmp(neuronTable.animalID,"190425") |...
%     strcmp(neuronTable.animalID,"190424") | strcmp(neuronTable.animalID,"190323") |...
%     strcmp(neuronTable.animalID,"190213"),:); %DR with light cabinet worry

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190426") |...
%     strcmp(neuronTable.animalID,"190424") | strcmp(neuronTable.animalID,"190323") |...
%     strcmp(neuronTable.animalID,"190213"),:); %DR with light cabinet worry

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190703") | strcmp(neuronTable.animalID,"190615") |...
%     strcmp(neuronTable.animalID,"190613"),:); %DR after we resolved the light cab worry

%% Loading the different groups

%ADULTS:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190116") | strcmp(neuronTable.animalID,"190117") |...
%     strcmp(neuronTable.animalID,"190404") | strcmp(neuronTable.animalID,"190405") |...
%     strcmp(neuronTable.animalID,"190419"),:); %The adults
% 
% neuronTable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"NR"),:);



%P14s:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190522") | strcmp(neuronTable.animalID,"190523"),:); %The devs
% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190523"),:); 
% 
% neuronTable = neuronTable(neuronTable.age < 16 & strcmp(neuronTable.condition,"NR"),:);




%DARK REARS:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190424") | strcmp(neuronTable.animalID,"190425") |...
%     strcmp(neuronTable.animalID,"190426") | strcmp(neuronTable.animalID,"190323")...
%     ,:); %The dark-rears

neuronTable = neuronTable(neuronTable.age > 29 & strcmp(neuronTable.condition,"DR"),:);



%B2KOs:

% neuronTable = neuronTable(strcmp(neuronTable.animalID,"190424") | strcmp(neuronTable.animalID,"190425") |...
%     strcmp(neuronTable.animalID,"190426") | strcmp(neuronTable.animalID,"190323")...
%     ,:); %The dark-rears

% neuronTable = b2koClass;



%% 

% Set a DSI or VARIANCE thresholds??
DSIthresh = 0; %Set to 0 for no Thresh
DSIsigThresh = 0.95; 
VARthresh = 10000; %Set to 10 for no Thresh; Set to 0.2 for reasonable thresh
dFoFthresh = 0;
neuronTable = neuronTable(neuronTable.DSI > DSIthresh & neuronTable.DSIsig > DSIsigThresh &...
    neuronTable.varSum < VARthresh & max(neuronTable.meanRespToBars,[],2) > dFoFthresh,:);





% % Do you just want to look at HB9s?
% neuronTable = neuronTable(strcmp(neuronTable.GFPid, "trhr") ,:);
% neuronTable = neuronTable(strcmp(neuronTable.GFPid, "drd4") ,:);
% neuronTable = neuronTable(strcmp(neuronTable.GFPid, "hb9") ,:);
% neuronTable = neuronTable(strcmp(neuronTable.cellID,"ON-OFF"),:);
neuronTable = neuronTable(strcmp(neuronTable.cellID,"ON"),:);



%% Retina length split in 4   - polar plot with single lines
VS_axis = 0.5;
% numBins = 24; %original bin size chosen by Alex
numBins = 36; %to match Remi

retinaFigBinnedInSpace = figure;
set(gcf, 'Position', [100 650 1600 150]);

% tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);
tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal")| strcmp(neuronTable.location,"ventroTemporal"),:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(1,10,1)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(1,10,2:3)
plot(mean(tempSil));
subplot(1,10,4:5)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal"),:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(1,10,6)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(1,10,7:8)
plot(mean(tempSil));
subplot(1,10,9:10)
histPlot(tempTable.prefDirCorr,numBins)

set(gcf, 'Renderer', 'painters');


retinaFigDiv4 = figure; 

set(gcf, 'Position', [100 50 1600 600]);
% sgtitle('placeHolder')

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON<2000/4,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,1)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,2:3)
plot(mean(tempSil));
subplot(4,10,4:5)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>2000/4 & neuronTable.distanceFromON<2*2000/4 ,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,11)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,12:13)
plot(mean(tempSil));
subplot(4,10,14:15)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>2*2000/4 & neuronTable.distanceFromON<3*2000/4 ,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,21)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,22:23)
plot(mean(tempSil));
subplot(4,10,24:25)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>3*2000/4,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,31)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,32:33)
plot(mean(tempSil));
subplot(4,10,34:35)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON<2000/4,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,6)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,7:8)
plot(mean(tempSil));
subplot(4,10,9:10)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>2000/4 & neuronTable.distanceFromON<2*2000/4 ,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,16)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,17:18)
plot(mean(tempSil));
subplot(4,10,19:20)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>2*2000/4 & neuronTable.distanceFromON<3*2000/4 ,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,26)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,27:28)
plot(mean(tempSil));
subplot(4,10,29:30)
histPlot(tempTable.prefDirCorr,numBins)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>3*2000/4,:);
tempSil = silTest(tempTable.prefDirCorr);
subplot(4,10,36)
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)
subplot(4,10,37:38)
plot(mean(tempSil));
subplot(4,10,39:40)
histPlot(tempTable.prefDirCorr,numBins)

set(gcf, 'Renderer', 'painters');

%% Figure with axis handles (All at top, split in 4 bottom)

VS_axis = 0.5;


tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);


imageSize = [100 100 900 900];
xLength = 1/4; % 4 columns in plot
yLength = 1/5; % 5 rows in the plot

fig_all = figure;
set(gcf, 'Position', imageSize);

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0.8 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0.8 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal"),:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0.8 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0.8 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON<2000/4,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0.6 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0.6 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>2000/4 & neuronTable.distanceFromON<2*2000/4 ,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0.4 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0.4 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>2*2000/4 & neuronTable.distanceFromON<3*2000/4,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0.2 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0.2 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>3*2000/4,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON<2000/4,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0.6 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0.6 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>2000/4 & neuronTable.distanceFromON<2*2000/4 ,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0.4 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0.4 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>2*2000/4 & neuronTable.distanceFromON<3*2000/4,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0.2 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0.2 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>3*2000/4,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

set(gcf, 'Renderer', 'painters');

%% Figure with axis handles (All at top, split in 2 bottom)

VS_axis = 0.5;


tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal"),:);


imageSize = [100 100 900 900];
xLength = 1/4; % 4 columns in plot
yLength = 1/5; % 5 rows in the plot

fig_all = figure;
set(gcf, 'Position', imageSize);

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON<2000/3,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0.66 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0.66 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON<2000/3,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0.66 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0.66 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON<2*2000/3 & neuronTable.distanceFromON>2000/3,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0.33 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0.33 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroNasal") & neuronTable.distanceFromON>2000/3,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0 0 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.25 0 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON<2*2000/3 & neuronTable.distanceFromON>2000/3,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0.33 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0.33 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

tempTable = neuronTable(strcmp(neuronTable.location,"ventroTemporal") & neuronTable.distanceFromON>2000/3,:);
axis_histPol = axes('Units', 'normalized', 'Position',[0.5 0 xLength yLength]);
DSpolHistPlot(tempTable.prefDirCorr);
axis_linePol = axes('Units', 'normalized', 'Position',[0.75 0 xLength yLength]);
DSpolPlot(tempTable.prefDirCorr,tempTable.vecSum, VS_axis)

set(gcf, 'Renderer', 'painters');


%% Functions

function DSpolPlot(prefDir,vecSum, VS_axis)

    polarplot([prefDir prefDir]', [zeros(length(prefDir),1) vecSum]','r');
    rlim([0 VS_axis]);
    set(gca,'thetaticklabel',{[]})
    set(gca,'rticklabel',{[]})
    thetaticks([])
    rticks([])

end

function DSpolHistPlot(prefDir)
%     polarhistogram(prefDir,40)
    numBins = 40;
    polarhistogram(prefDir,'BinEdges',0:pi/numBins*2:2*pi)
%     set(gca,'thetaticklabel',{[]})
%     set(gca,'rticklabel',{[]})
    thetaticks([0 45 90 135 180 225 270 315])
%     rticks([])
%     rlim([0 70])

end


function histPlot(prefDir,numBins)

    histogram(prefDir,'BinEdges',0:pi/numBins*2:2*pi);
%     ylim([0 20])
%     xlim([0 2*pi])
    xticks(0:pi/4:2*pi)
    xticklabels({'0','45','90','135','180','225','270','315','360'})

end
