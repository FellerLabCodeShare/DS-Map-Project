load ventralAllClusteredWithLocs.mat
load b2koClass.mat

% Set a DSI or VARIANCE thresholds??
DSIthresh = 0; %Set to 0 for no Thresh
DSIsigThresh = 0.95; 
VARthresh = 10000; %Set to 10 for no Thresh; Set to 0.2 for reasonable thresh
dFoFthresh = 0;
b2koClass = b2koClass(b2koClass.DSI > DSIthresh & b2koClass.DSIsig > DSIsigThresh &...
    b2koClass.varSum < VARthresh & max(b2koClass.meanRespToBars,[],2) > dFoFthresh &...
    b2koClass.age > 29,:);

%% add exp grps to tables

EO.expGrp = repmat("EO",[height(EO),1]);
NR.expGrp = repmat("NR",[height(NR),1]);
DR.expGrp = repmat("DR",[height(DR),1]);
b2koClass.expGrp = repmat("b2koClass",[height(b2koClass),1]);

allTable = table;
allTable.DSI = [EO.DSI;NR.DSI;DR.DSI;b2koClass.DSI];
allTable.vecSum = [EO.vecSum;NR.vecSum;DR.vecSum;b2koClass.vecSum];
allTable.expID = [EO.expGrp;NR.expGrp;DR.expGrp;b2koClass.expGrp];
allTable.cellID = [EO.cellID; NR.cellID; DR.cellID; b2koClass.cellID];
allTable.dsID = [EO.idxDS;NR.idxDS;DR.idxDS;b2koClass.idxDS];
allTable.orientation = repmat("LOL",[height(allTable),1]);


% onOffTable = allTable(strcmp(allTable.cellID,"ON-OFF"),:);
% onTable = allTable(strcmp(allTable.cellID,"ON"),:);

%% get index for different subgroups

ind_allHor = allTable.dsID == 1 | allTable.dsID == 3;
ind_allVer = allTable.dsID == 2 | allTable.dsID == 4;

ind_onOffHor = allTable.dsID == 1 | allTable.dsID == 3 & strcmp(allTable.cellID,"ON-OFF");
ind_onOffVer = allTable.dsID == 2 | allTable.dsID == 4 & strcmp(allTable.cellID,"ON-OFF");

ind_onHor = allTable.dsID == 1 | allTable.dsID == 3 & strcmp(allTable.cellID,"ON");
ind_onVer = allTable.dsID == 2 | allTable.dsID == 4 & strcmp(allTable.cellID,"ON");

allTable.orientation(ind_allHor) = "Hor";
allTable.orientation(ind_allVer) = "Ver";


%% plot DSIs

%all DSGCs combined
figure
boxplot(allTable.DSI(ind_allHor), allTable.expID(ind_allHor),'Notch','on');
ylim([0 1]);
figure
boxplot(allTable.DSI(ind_allVer), allTable.expID(ind_allVer),'Notch','on');
ylim([0 1]);

% tempTable = allTable(ind_allHor,:);
tempTable = allTable(ind_allVer,:);
figure('Name','AllHor'), hold
ecdf(tempTable.DSI(strcmp(tempTable.expID,"NR")));
ecdf(tempTable.DSI(strcmp(tempTable.expID,"EO")));
ecdf(tempTable.DSI(strcmp(tempTable.expID,"DR")));
ecdf(tempTable.DSI(strcmp(tempTable.expID,"b2koClass")));
xlim([0 1])

%on-offs only
figure
boxplot(allTable.DSI(ind_onOffHor), allTable.expID(ind_onOffHor),'Notch','on');
ylim([0 1]);
figure
boxplot(allTable.DSI(ind_onOffVer), allTable.expID(ind_onOffVer),'Notch','on');
ylim([0 1]);

%ons only
figure
boxplot(allTable.DSI(ind_onHor), allTable.expID(ind_onHor),'Notch','on');
ylim([0 1]);
figure
boxplot(allTable.DSI(ind_onVer), allTable.expID(ind_onVer),'Notch','on');
ylim([0 1]);


% Plot DSIs as violin plots
figure, violin({allTable.DSI(strcmp(allTable.expID,"EO") & strcmp(allTable.orientation,"Hor")),allTable.DSI(strcmp(allTable.expID,"NR") & strcmp(allTable.orientation,"Hor")),...
    allTable.DSI(strcmp(allTable.expID,"DR") & strcmp(allTable.orientation,"Hor")),allTable.DSI(strcmp(allTable.expID,"b2koClass") & strcmp(allTable.orientation,"Hor"))})
ylim([0 1])
figure, violin({allTable.DSI(strcmp(allTable.expID,"EO") & strcmp(allTable.orientation,"Ver")),allTable.DSI(strcmp(allTable.expID,"NR") & strcmp(allTable.orientation,"Ver")),...
    allTable.DSI(strcmp(allTable.expID,"DR") & strcmp(allTable.orientation,"Ver")),allTable.DSI(strcmp(allTable.expID,"b2koClass") & strcmp(allTable.orientation,"Ver"))})
ylim([0 1])




%% plot vecSums
% 
% %all DSGCs combined
% figure
% boxplot(allTable.vecSum(ind_allHor), allTable.expID(ind_allHor),'Notch','on');
% ylim([0 1]);
% figure
% boxplot(allTable.vecSum(ind_allVer), allTable.expID(ind_allVer),'Notch','on');
% ylim([0 1]);
% 
% %on-offs only
% figure
% boxplot(allTable.vecSum(ind_onOffHor), allTable.expID(ind_onOffHor),'Notch','on');
% ylim([0 1]);
% figure
% boxplot(allTable.vecSum(ind_onOffVer), allTable.expID(ind_onOffVer),'Notch','on');
% ylim([0 1]);
% 
% %ons only
% figure
% boxplot(allTable.vecSum(ind_onHor), allTable.expID(ind_onHor),'Notch','on');
% ylim([0 1]);
% figure
% boxplot(allTable.vecSum(ind_onVer), allTable.expID(ind_onVer),'Notch','on');
% ylim([0 1]);


% % Plot DSIs as violin plots
% figure, violin({allTable.vecSum(strcmp(allTable.expID,"EO") & strcmp(allTable.orientation,"Hor")),allTable.vecSum(strcmp(allTable.expID,"NR") & strcmp(allTable.orientation,"Hor")),...
%     allTable.vecSum(strcmp(allTable.expID,"DR") & strcmp(allTable.orientation,"Hor")),allTable.vecSum(strcmp(allTable.expID,"b2koClass") & strcmp(allTable.orientation,"Hor"))})
% ylim([0 1])
% figure, violin({allTable.vecSum(strcmp(allTable.expID,"EO") & strcmp(allTable.orientation,"Ver")),allTable.vecSum(strcmp(allTable.expID,"NR") & strcmp(allTable.orientation,"Ver")),...
%     allTable.vecSum(strcmp(allTable.expID,"DR") & strcmp(allTable.orientation,"Ver")),allTable.vecSum(strcmp(allTable.expID,"b2koClass") & strcmp(allTable.orientation,"Ver"))})
% ylim([0 1])



% %% now stats
% 
% %%%% here is the stats for runing one way anovas on 4 exp grp of ON or
% %%%% ON-OFF within each cell subtype
% 
% %comment the one you want: ON-0FF or ON
% % tempTable = allTable(strcmp(allTable.cellID,"ON-OFF"),:);
% tempTable = allTable(strcmp(allTable.cellID,"ON"),:);
% 
% % tempTable = tempTable(tempTable.dsID == 1,:);
% % tempTable = tempTable(tempTable.dsID == 2,:);
% % tempTable = tempTable(tempTable.dsID == 3,:);
% tempTable = tempTable(tempTable.dsID == 4,:);
% 
% %Now doing posthoc to compare exp grps
% 
% %1 way anova on single cell types
% [p,table,stats] = anova1(tempTable.DSI,tempTable.expID);
% results = multcompare(stats)


% %%%% here is the stats for runing one way anovas on 4 subtypes of ON or
% %%%% ON-OFF within each exp group
% 
% %comment the one you want: ON-0FF or ON
% % tempTable = allTable(strcmp(allTable.cellID,"ON-OFF"),:);
% tempTable = allTable(strcmp(allTable.cellID,"ON"),:);
% 
% tempTable = tempTable(strcmp(tempTable.expID,"EO"),:);
% % tempTable = tempTable(strcmp(tempTable.expID,"NR"),:);
% % tempTable = tempTable(strcmp(tempTable.expID,"DR"),:);
% % tempTable = tempTable(strcmp(tempTable.expID,"b2koClass"),:);
% 
% 
% %Now doing posthoc to compare cell types to each other
% 
% %1 way anova on single cell types
% [p,table,stats] = anova1(tempTable.vecSum,tempTable.dsID);
% results = multcompare(stats)


% % % % % going to reformart data so I can feed into a 2 way anova. One factor will
% % % % % be horizontal vs vertical, the other NR or b2ko
% 
% % % % %Doing separate anovas for all ds, On off DS, and ON ds
% % testTable = allTable; %For both ON OFF and ON
% % testTable = allTable(strcmp(allTable.cellID,"ON-OFF"),:);
% testTable = allTable(strcmp(allTable.cellID,"ON"),:);
% 
% testTable = testTable(strcmp(testTable.expID,"NR") | strcmp(testTable.expID,"DR"),:);
% 
% % % % % 2 factor anova (separating ON-OFF and ON).
% [p,table,stats] = anovan(testTable.vecSum,{testTable.expID testTable.orientation},'model','full',...
%     'varnames',{'expID','orientation'});
% 
% results = multcompare(stats,'Dimension',[1 2])
% 
% 
% 
% % 3 factor anova (Add ON-OFF and ON as 3rd level)
% % [p,table,stats] = anovan(testTable.vecSum,{testTable.expID testTable.orientation testTable.cellID},'model','full',...
% %     'varnames',{'expID','orientation','cellID'});
% % 
% % results = multcompare(stats,'Dimension',[1 3])











