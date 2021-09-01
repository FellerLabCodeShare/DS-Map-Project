
load ventralAllClusteredWithLocs.mat
load b2koClass.mat

% Set a DSI or VARIANCE thresholds??
DSIthresh = 0; %Set to 0 for no Thresh
DSIsigThresh = 0.95; 
VARthresh = 10000; %Set to 10 for no Thresh; Set to 0.2 for reasonable thresh
dFoFthresh = 0;

b2koEO = b2koClass(strcmp(b2koClass.animalID,"201007") | strcmp(b2koClass.animalID,"201028") |...
    strcmp(b2koClass.animalID,"201223") | strcmp(b2koClass.animalID,"210111") |...
    strcmp(b2koClass.animalID,"210112"),:); %EO B2KOS
b2koEO = b2koEO(b2koEO.DSI > DSIthresh & b2koEO.DSIsig > DSIsigThresh &...
    b2koEO.varSum < VARthresh & max(b2koEO.meanRespToBars,[],2) > dFoFthresh &...
    b2koEO.age <20,:);


b2koClass = b2koClass(b2koClass.DSI > DSIthresh & b2koClass.DSIsig > DSIsigThresh &...
    b2koClass.varSum < VARthresh & max(b2koClass.meanRespToBars,[],2) > dFoFthresh &...
    b2koClass.age >20,:);




% NR = NR(strcmp(NR.cellID,"ON-OFF"),:);
% DR = DR(strcmp(DR.cellID,"ON-OFF"),:);
% EO = EO(strcmp(EO.cellID,"ON-OFF"),:);

% %% ventroNasal vs ventroTemporal
% 
% 
% NR_VN = NR(strcmp(NR.location,"ventroNasal"),:);
% NR_VT = NR(strcmp(NR.location,"ventroTemporal"),:);
% 
% plotProps(NR_VN)
% plotProps(NR_VT)
% 
% plotDSI(NR_VN)
% plotDSI(NR_VT)
% 
% 
% %% Central vs Peripheral
% 
% NR_central = NR(NR.distanceFromON < 1000,:);
% NR_peripheral = NR(NR.distanceFromON > 999,:);
% 
% plotProps(NR_central)
% plotProps(NR_peripheral)
% 
% plotDSI(NR_central)
% plotDSI(NR_peripheral)

%% Proportions

%First plot
[EOpropRes, EOcellID, EOpropResHorVer, EOcellIDhorVer] = plotProps(EO);
[NRpropRes, NRcellID, NRpropResHorVer, NRcellIDhorVer] = plotProps(NR);
[DRpropRes, DRcellID, DRpropResHorVer, DRcellIDhorVer] = plotProps(DR);
[b2kopropRes, b2kocellID, b2kopropResHorVer, b2kocellIDhorVer] = plotProps(b2koClass);
[b2koEOpropRes, b2koEOcellID, b2koEOpropResHorVer, b2koEOcellIDhorVer] = plotProps(b2koEO);

% %GFP stats
% [dummy, dummy, dummy, dummy,EOpropsGFP,EOgfpID] = plotProps(EO);
% [dummy, dummy, dummy, dummy,NRpropsGFP,NRgfpID] = plotProps(NR);
% [dummy, dummy, dummy, dummy,DRpropsGFP,DRgfpID] = plotProps(DR);
% 
% gfpTable = table;
% gfpTable.props = [EOpropsGFP;NRpropsGFP;DRpropsGFP];
% gfpTable.gfpID = [EOgfpID;NRgfpID;DRgfpID];
% gfpTable.expID = [repmat("EO",size(EOgfpID));repmat("NR",size(NRgfpID));repmat("DR",size(DRgfpID))];
% gfpTable = gfpTable(strcmp(gfpTable.gfpID,"drd4"),:);
% [p,tbl,stats] = anova1(gfpTable.props,gfpTable.expID);
% c = multcompare(stats);

%Reformat proportions that were calculated per FOV

propsAll = [EOpropRes;NRpropRes;DRpropRes;b2kopropRes;b2koEOpropRes];
expID = [repmat("EO",size(EOcellID));repmat("NR",size(NRcellID));repmat("DR",size(DRcellID));repmat("b2ko",size(b2kocellID));repmat("b2koEO",size(b2koEOcellID))];
cellIDall = [EOcellID;NRcellID;DRcellID;b2kocellID;b2koEOcellID];
%strings in cellIDall: "ooT" "ooD" "ooN" "ooV" "oT" "oD" "oN" "oV"

propsTableSubtypes = table;
propsTableSubtypes.props = propsAll;
propsTableSubtypes.expID = expID;
propsTableSubtypes.stringID = cellIDall;
propsTableSubtypes.onOffID(strcmp(propsTableSubtypes.stringID,'ooT') | strcmp(propsTableSubtypes.stringID,'ooN')| strcmp(propsTableSubtypes.stringID,'ooV')| strcmp(propsTableSubtypes.stringID,'ooD')) = "ON-OFF";
propsTableSubtypes.onOffID(strcmp(propsTableSubtypes.stringID,'oT') | strcmp(propsTableSubtypes.stringID,'oN')| strcmp(propsTableSubtypes.stringID,'oV')| strcmp(propsTableSubtypes.stringID,'oD')) = "ON";
% propsTableSubtypes.orientation(strcmp(propsTable.stringID,'ooHor') | strcmp(propsTable.stringID,'oHor')) = "Hor";
% propsTableSubtypes.orientation(strcmp(propsTable.stringID,'ooVer') | strcmp(propsTable.stringID,'oVer')) = "Ver";

propsAllHorVer = [EOpropResHorVer;NRpropResHorVer;DRpropResHorVer;b2kopropResHorVer;b2koEOpropResHorVer];
expIDHorVer = [repmat("EO",size(EOcellIDhorVer));repmat("NR",size(NRcellIDhorVer));repmat("DR",size(DRcellIDhorVer));repmat("b2ko",size(b2kocellIDhorVer));repmat("b2koEO",size(b2koEOcellIDhorVer))];
cellIDhorVer = [EOcellIDhorVer;NRcellIDhorVer;DRcellIDhorVer;b2kocellIDhorVer;b2koEOcellIDhorVer];
%strings in cellIDhorVer: "allHor" "allVer" "ooHor" "ooVer" "oHor" "oVer"

propsTable = table;
propsTable.props = propsAllHorVer;
propsTable.expID = expIDHorVer;
propsTable.stringID = cellIDhorVer;
propsTable.onOffID(strcmp(propsTable.stringID,'ooHor') | strcmp(propsTable.stringID,'ooVer')) = "ON-OFF";
propsTable.onOffID(strcmp(propsTable.stringID,'oHor') | strcmp(propsTable.stringID,'oVer')) = "ON";
propsTable.orientation(strcmp(propsTable.stringID,'ooHor') | strcmp(propsTable.stringID,'oHor')) = "Hor";
propsTable.orientation(strcmp(propsTable.stringID,'ooVer') | strcmp(propsTable.stringID,'oVer')) = "Ver";


%% CDF plots

% subtype = 'ooN';
% tempTable = propsTableSubtypes(strcmp(propsTableSubtypes.stringID,subtype),:); %Per subtype
% tempTable = propsTable(strcmp(propsTable.stringID,subtype),:); %all hor/ver
% figure, hold
% % ecdf(tempTable.props(strcmp(tempTable.expID,'EO')));
% ecdf(tempTable.props(strcmp(tempTable.expID,'NR')));
% % ecdf(tempTable.props(strcmp(tempTable.expID,'DR')));
% ecdf(tempTable.props(strcmp(tempTable.expID,'b2ko')))
% ecdf(tempTable.props(strcmp(tempTable.expID,'b2koEO')))
% xlim([0 0.6])

% tempTable(strcmp(tempTable.expID,"DR"),:) = [];
% tempTable(strcmp(tempTable.expID,"b2ko"),:) = [];
% tempTable(strcmp(tempTable.expID,"b2koEO"),:) = [];

% [p,tbl,stats] = anova1(tempTable.props,tempTable.expID);
% c = multcompare(stats);


b2table = propsTableSubtypes(strcmp(propsTableSubtypes.expID,"b2koEO"),:);
b2table = b2table(strcmp(b2table.onOffID,"ON"),:);
[p,tbl,stats] = anova1(b2table.props,b2table.stringID);
c = multcompare(stats);


% figure, hold
% violin({tempTable.props(strcmp(tempTable.expID,'EO')),tempTable.props(strcmp(tempTable.expID,'NR')),tempTable.props(strcmp(tempTable.expID,'DR')),tempTable.props(strcmp(tempTable.expID,'b2ko'))})
% ylim([0 0.6])
% 
% figure, hold
% plotSpread({tempTable.props(strcmp(tempTable.expID,'EO')),tempTable.props(strcmp(tempTable.expID,'NR')),tempTable.props(strcmp(tempTable.expID,'DR')),tempTable.props(strcmp(tempTable.expID,'b2ko'))})
% ylim([0 0.6])


% figure, boxplot(propsAllHorVer(strcmp(cellIDhorVer,"allHor")),expIDHorVer(strcmp(cellIDhorVer,"allHor")),'Notch','on');
% axis([0 5 0 1])
% figure, boxplot(propsAllHorVer(strcmp(cellIDhorVer,"allVer")),expIDHorVer(strcmp(cellIDhorVer,"allVer")),'Notch','on');
% axis([0 5 0 1])
% figure, boxplot(propsAllHorVer(strcmp(cellIDhorVer,"ooHor")),expIDHorVer(strcmp(cellIDhorVer,"ooHor")),'Notch','on');
% axis([0 5 0 1])
% figure, boxplot(propsAllHorVer(strcmp(cellIDhorVer,"ooVer")),expIDHorVer(strcmp(cellIDhorVer,"ooVer")),'Notch','on');
% axis([0 5 0 1])
% figure, boxplot(propsAllHorVer(strcmp(cellIDhorVer,"oHor")),expIDHorVer(strcmp(cellIDhorVer,"oHor")),'Notch','on');
% axis([0 5 0 1])
% figure, boxplot(propsAllHorVer(strcmp(cellIDhorVer,"oVer")),expIDHorVer(strcmp(cellIDhorVer,"oVer")),'Notch','on');
% axis([0 5 0 1])

% figure, violin({propsTable.props(strcmp(propsTable.expID,"EO") & strcmp(propsTable.stringID,"allHor")),propsTable.props(strcmp(propsTable.expID,"NR") & strcmp(propsTable.stringID,"allHor")),...
%     propsTable.props(strcmp(propsTable.expID,"DR") & strcmp(propsTable.stringID,"allHor")),propsTable.props(strcmp(propsTable.expID,"b2ko") & strcmp(propsTable.stringID,"allHor"))})
% ylim([0 1])
% figure, violin({propsTable.props(strcmp(propsTable.expID,"EO") & strcmp(propsTable.stringID,"allVer")),propsTable.props(strcmp(propsTable.expID,"NR") & strcmp(propsTable.stringID,"allVer")),...
%     propsTable.props(strcmp(propsTable.expID,"DR") & strcmp(propsTable.stringID,"allVer")),propsTable.props(strcmp(propsTable.expID,"b2ko") & strcmp(propsTable.stringID,"allVer"))})
% ylim([0 1])



% %% stats for subtypes
% %comment the one you want: ON-0FF or ON
% % tempTable = propsTableSubtypes(strcmp(propsTableSubtypes.onOffID,"ON-OFF"),:);
% tempTable = propsTableSubtypes(strcmp(propsTableSubtypes.onOffID,"ON"),:);
% % 
% % tempTable = tempTable(strcmp(tempTable.expID,"EO"),:);
% % tempTable = tempTable(strcmp(tempTable.expID,"NR"),:);
% % tempTable = tempTable(strcmp(tempTable.expID,"DR"),:);
% tempTable = tempTable(strcmp(tempTable.expID,"b2ko"),:);
% 
% 
% %Now doing posthoc comparisons within each cell type
% 
% %1 way anova on single cell types
% [p,table,stats] = anova1(tempTable.props,tempTable.stringID);
% results = multcompare(stats)


%% stats for HOR VS VER
% %comment the one you want: ON-0FF or ON
% tempTable = propsTable(strcmp(propsTable.onOffID,"ON-OFF"),:);
% % tempTable = propsTable(strcmp(propsTable.onOffID,"ON"),:);
% 
% tempTable = tempTable(strcmp(tempTable.expID,"NR") | strcmp(tempTable.expID,"b2ko"),:);
% 
% tempTable(strcmp(tempTable.stringID,"allHor") | strcmp(tempTable.stringID,"allVer"),:) = [];

% % % % % % % [p,table,stats] = anovan(tempTable.props,{tempTable.expID tempTable.orientation tempTable.onOffID},'model','full',...
% % % % % % %     'varnames',{'expID','orientation','onOffID'});



% [p,table,stats] = anovan(tempTable.props,{tempTable.expID tempTable.orientation},'model','full',...
%     'varnames',{'expID','orientation'});
% 
% results = multcompare(stats,'Dimension',[1 2])

% 

% % % % % % 
% % % % % % % %2 factor anova on whole dataset
% % % % % % % [p,table,stats] = anovan(propsAll,{cellIDall expID},'model','interaction',...
% % % % % % %     'varnames',{'cellIDall','expID'});

%Results show main effect of cell type but not of experimental condition.
%Effect in the interaction of the two also present, so...



% figure, boxplot(propsAll(ind_ALLhori),expID(ind_ALLhori),'Notch','on');
% axis([0 5 0 0.6])




%Now doing posthoc comparisons within each cell type

% %1 way anova on single cell types
% ooTind = strcmp(cellIDall,'ooT');
% [p,table,stats] = anova1(propsAll(ooTind),expID(ooTind));








%% Tuning

EO.expGrp = repmat("EO",[height(EO),1]);
NR.expGrp = repmat("NR",[height(NR),1]);
DR.expGrp = repmat("DR",[height(DR),1]);
b2koClass.expGrp = repmat("b2ko",[height(b2koClass),1]);
b2koEO.expGrp = repmat("b2koEO",[height(b2koEO),1]);

tuningTable = table;
tuningTable.DSI = [EO.DSI;NR.DSI;DR.DSI;b2koClass.DSI;b2koEO.DSI];
tuningTable.vecSum = [EO.vecSum;NR.vecSum;DR.vecSum;b2koClass.vecSum;b2koEO.vecSum];
tuningTable.expID = [EO.expGrp;NR.expGrp;DR.expGrp;b2koClass.expGrp;b2koEO.expGrp];
tuningTable.cellID = [EO.cellID; NR.cellID; DR.cellID; b2koClass.cellID;b2koEO.cellID];
tuningTable.dsID = [EO.idxDS;NR.idxDS;DR.idxDS;b2koClass.idxDS;b2koEO.idxDS];
tuningTable.orientation = repmat("LOL",[height(tuningTable),1]);
tuningTable.gfpID = [EO.GFPid;NR.GFPid;DR.GFPid;b2koClass.GFPid;b2koEO.GFPid];


ind_allHor = tuningTable.dsID == 1 | tuningTable.dsID == 3;
ind_allVer = tuningTable.dsID == 2 | tuningTable.dsID == 4;

ind_onOffHor = tuningTable.dsID == 1 | tuningTable.dsID == 3 & strcmp(tuningTable.cellID,"ON-OFF");
ind_onOffVer = tuningTable.dsID == 2 | tuningTable.dsID == 4 & strcmp(tuningTable.cellID,"ON-OFF");

ind_onHor = tuningTable.dsID == 1 | tuningTable.dsID == 3 & strcmp(tuningTable.cellID,"ON");
ind_onVer = tuningTable.dsID == 2 | tuningTable.dsID == 4 & strcmp(tuningTable.cellID,"ON");

tuningTable.orientation(ind_allHor) = "Hor";
tuningTable.orientation(ind_allVer) = "Ver";


% Stats for DSI or VS

%For Figure 1
table4stats = tuningTable;
table4stats = table4stats(strcmp(table4stats.cellID,"ON"),:);

% table4stats = table4stats(table4stats.dsID == 1,:);
% table4stats = table4stats(table4stats.dsID == 2 | table4stats.dsID == 4,:);
% table4stats = table4stats(strcmp(table4stats.gfpID,"drd4"),:);
% table4stats(strcmp(table4stats.expID,"DR"),:) = [];
% table4stats(strcmp(table4stats.expID,"b2ko"),:) = [];
% table4stats(strcmp(table4stats.expID,"b2koEO"),:) = [];
% [p,tbl,stats] = anova1(table4stats.vecSum,table4stats.expID);
% c = multcompare(stats);

tempTable = table4stats(strcmp(table4stats.expID,"b2koEO"),:);
% [p,tbl,stats] = anova1(tempTable.vecSum,tempTable.dsID);
% c = multcompare(stats);
% 
% SEM = std(tempTable.DSI)/sqrt(length(tempTable.DSI));       % Standard Error
% ts = tinv([0.025  0.975],length(tempTable.DSI)-1);          % T-Score
% CI = mean(tempTable.DSI) + ts*SEM;                          % Confidence Intervals



% %For hor and ver (figure 2)
% table4stats = tuningTable;
% table4stats = table4stats(strcmp(table4stats.orientation,"Hor"),:);
% table4stats(strcmp(table4stats.expID,"DR"),:) = [];
% [p,tbl,stats] = anova1(table4stats.DSI,table4stats.expID);
% c = multcompare(stats);










% % 
% plotDSI(EO)
% plotDSI(NR)
% plotDSI(DR)
% plotDSI(b2koClass)
% plotDSI(b2koEO)

% 
% % 
% plotVS(EO)
% plotVS(NR)
% plotVS(DR)
% plotVS(b2koClass)
% plotVS(b2koEO)
% % 
% % 





function [propResColumnar, cellID, propResultsConcHorVerColumnar, cellIDconcHorVer] = plotProps(inputTable)
% function [propResColumnar, cellID, propResultsConcHorVerColumnar, cellIDconcHorVer, propResGFPColumnar, cellIDGFP] = plotProps(inputTable)

numFOVs = unique(inputTable.fileName);
hb9FOVs = unique(inputTable.fileName(strcmp(inputTable.GFPretina,'hb9')));
drd4FOVs = unique(inputTable.fileName(strcmp(inputTable.GFPretina,'drd4')));

ooTable = inputTable(strcmp(inputTable.cellID,"ON-OFF"),:);
oTable = inputTable(strcmp(inputTable.cellID,"ON"),:);

all_hor = [];
all_ver = [];
oo_hor = [];
oo_ver = [];
o_hor = [];
o_ver = [];

oo_numT = [];
oo_numD = [];
oo_numN = [];
oo_numV = [];
o_numT = [];
o_numD = [];
o_numN = [];
o_numV = [];
hb9_num = [];
drd4_num = [];

    for i = 1:length(numFOVs)
        oo_tempTable = ooTable(strcmp(ooTable.fileName,numFOVs(i)),:);
        o_tempTable = oTable(strcmp(oTable.fileName,numFOVs(i)),:);

        oo_numNeurons = length(oo_tempTable.idxDS);
        o_numNeurons = length(o_tempTable.idxDS);
        totalNeurons = oo_numNeurons+o_numNeurons;
        
        if sum(strcmp(numFOVs(i),hb9FOVs)) == 1
            hb9_num = [hb9_num; sum(strcmp(oo_tempTable.GFPid,'hb9'))/totalNeurons];
        else 
            hb9_num = [hb9_num; NaN];
        end
        
        if sum(strcmp(numFOVs(i),drd4FOVs)) == 1
            drd4_num = [drd4_num; sum(strcmp(oo_tempTable.GFPid,'drd4'))/totalNeurons];
        else 
            drd4_num = [drd4_num; NaN];
        end
        
        all_hor = [all_hor; (length(oo_tempTable.idxDS(oo_tempTable.idxDS == 1 | oo_tempTable.idxDS == 3)) + length(o_tempTable.idxDS(o_tempTable.idxDS == 1 | o_tempTable.idxDS == 3)))/totalNeurons];
        all_ver = [all_ver; (length(oo_tempTable.idxDS(oo_tempTable.idxDS == 2 | oo_tempTable.idxDS == 4)) + length(o_tempTable.idxDS(o_tempTable.idxDS == 2 | o_tempTable.idxDS == 4)))/totalNeurons];

        oo_hor = [oo_hor;(length(oo_tempTable.idxDS(oo_tempTable.idxDS == 1 | oo_tempTable.idxDS == 3)))/totalNeurons];
        oo_ver = [oo_ver;(length(oo_tempTable.idxDS(oo_tempTable.idxDS == 2 | oo_tempTable.idxDS == 4)))/totalNeurons];
        o_hor = [o_hor;(length(o_tempTable.idxDS(o_tempTable.idxDS == 1 | o_tempTable.idxDS == 3)))/totalNeurons];
        o_ver = [o_ver;(length(o_tempTable.idxDS(o_tempTable.idxDS == 2 | o_tempTable.idxDS == 4)))/totalNeurons];
        
        oo_numT = [oo_numT; sum(oo_tempTable.idxDS == 1)/totalNeurons];
        oo_numD = [oo_numD; sum(oo_tempTable.idxDS == 2)/totalNeurons];
        oo_numN = [oo_numN; sum(oo_tempTable.idxDS == 3)/totalNeurons];
        oo_numV = [oo_numV; sum(oo_tempTable.idxDS == 4)/totalNeurons];

        o_numT = [o_numT; sum(o_tempTable.idxDS == 1)/totalNeurons];
        o_numD = [o_numD; sum(o_tempTable.idxDS == 2)/totalNeurons];
        o_numN = [o_numN; sum(o_tempTable.idxDS == 3)/totalNeurons];
        o_numV = [o_numV; sum(o_tempTable.idxDS == 4)/totalNeurons];

    end

% figure, plotSpread({oo_numT, oo_numD, oo_numN, oo_numV, o_numT, o_numD, o_numN, o_numV},'showMM',2)
figure, boxplot([oo_numT;oo_numD;oo_numN;oo_numV;o_numT;o_numD;o_numN;o_numV;hb9_num;drd4_num],...
    [ones(length(oo_numT),1);ones(length(oo_numD),1)*2;ones(length(oo_numN),1)*3;ones(length(oo_numD),1)*4;...
    ones(length(o_numT),1)*5;ones(length(o_numD),1)*6;ones(length(o_numN),1)*7;ones(length(o_numV),1)*8;...
    ones(length(hb9_num),1)*9;ones(length(drd4_num),1)*10],'Notch','on')
axis([0 11 0 0.6])

figure, violin({oo_numT,oo_numD,oo_numN,oo_numV,o_numT,o_numD,o_numN,o_numV})
ylim([0 0.6])



propResults = [oo_numT,oo_numD,oo_numN,oo_numV,o_numT,o_numD,o_numN,o_numV];
[row,col] = size(propResults);
propResColumnar = reshape(propResults,[row*col,1]);
cellID = [repmat("ooT",[row,1]);repmat("ooD",[row,1]);repmat("ooN",[row,1]);repmat("ooV",[row,1]);...
    repmat("oT",[row,1]);repmat("oD",[row,1]);repmat("oN",[row,1]);repmat("oV",[row,1])];

propResultsConcHorVer = [all_hor, all_ver, oo_hor, oo_ver, o_hor, o_ver];
[row,col] = size(propResultsConcHorVer);
propResultsConcHorVerColumnar = reshape(propResultsConcHorVer,[row*col,1]);
cellIDconcHorVer = [repmat("allHor",[row,1]);repmat("allVer",[row,1]);repmat("ooHor",[row,1]);repmat("ooVer",[row,1]);...
    repmat("oHor",[row,1]);repmat("oVer",[row,1])];

propResultsGFPcells = [hb9_num,drd4_num];
[row,col] = size(propResultsGFPcells);
propResGFPColumnar = reshape(propResultsGFPcells,[row*col,1]);
cellIDGFP = [repmat("hb9",[row,1]);repmat("drd4",[row,1])];

% figure('Name','ON-OFFs props'), hold
% ecdf(oo_numN);
% ecdf(oo_numT);
% ecdf(oo_numD);
% ecdf(oo_numV);
% xlim([0 0.6])
% 
% figure('Name','ONs props'), hold
% ecdf(o_numN);
% ecdf(o_numT);
% ecdf(o_numD);
% ecdf(o_numV);
% xlim([0 0.6])

% figure('Name','hb9 and drd4 props'), hold
% ecdf(hb9_num);
% ecdf(drd4_num);
% xlim([0 0.6])

end

function plotDSI(inputTable)

ooTable = inputTable(strcmp(inputTable.cellID,"ON-OFF"),:);
oTable = inputTable(strcmp(inputTable.cellID,"ON"),:);

oo_T = ooTable.DSI(ooTable.idxDS == 1);
oo_D = ooTable.DSI(ooTable.idxDS == 2);
oo_N = ooTable.DSI(ooTable.idxDS == 3);
oo_V = ooTable.DSI(ooTable.idxDS == 4);
o_T = oTable.DSI(oTable.idxDS == 1);
o_D = oTable.DSI(oTable.idxDS == 2);
o_N = oTable.DSI(oTable.idxDS == 3);
o_V = oTable.DSI(oTable.idxDS == 4);

hb9s = ooTable.DSI(strcmp(ooTable.GFPid,'hb9'));
drd4s = ooTable.DSI(strcmp(ooTable.GFPid,'drd4'));

% figure, plotSpread({oo_numT, oo_numD, oo_numN, oo_numV, o_numT, o_numD, o_numN, o_numV},'showMM',2)
figure, boxplot([oo_T;oo_D;oo_N;oo_V;o_T;o_D;o_N;o_V;hb9s;drd4s],...
    [ones(length(oo_T),1);ones(length(oo_D),1)*2;ones(length(oo_N),1)*3;ones(length(oo_V),1)*4;...
    ones(length(o_T),1)*5;ones(length(o_D),1)*6;ones(length(o_N),1)*7;ones(length(o_V),1)*8;...
    ones(length(hb9s),1)*9; ones(length(drd4s),1)*10],'Notch','on')
axis([0 11 0 1])

figure, violin({oo_N,oo_T,oo_D,oo_V,o_N,o_T,o_D,o_V})
ylim([0 1])

% figure('Name','Horizontal v Vertical DS'), hold
% ecdf([oo_N;oo_T;o_N;o_T])
% ecdf([oo_V;oo_D;o_V;o_D])
% xlim([0 1])

% figure('Name','Hb9 vs Drd4 DSI'), hold
% ecdf(hb9s)
% ecdf(drd4s)
% xlim([0 1])

figure('Name','ON-OFFs'), hold
ecdf(oo_N);
ecdf(oo_T);
ecdf(oo_D);
ecdf(oo_V);
xlim([0 1])

figure('Name','ONs'), hold
ecdf(o_N);
ecdf(o_T);
ecdf(o_D);
ecdf(o_V);
xlim([0 1])

end

function plotVS(inputTable)

ooTable = inputTable(strcmp(inputTable.cellID,"ON-OFF"),:);
oTable = inputTable(strcmp(inputTable.cellID,"ON"),:);

oo_T = ooTable.vecSum(ooTable.idxDS == 1);
oo_D = ooTable.vecSum(ooTable.idxDS == 2);
oo_N = ooTable.vecSum(ooTable.idxDS == 3);
oo_V = ooTable.vecSum(ooTable.idxDS == 4);
o_T = oTable.vecSum(oTable.idxDS == 1);
o_D = oTable.vecSum(oTable.idxDS == 2);
o_N = oTable.vecSum(oTable.idxDS == 3);
o_V = oTable.vecSum(oTable.idxDS == 4);

hb9s = ooTable.vecSum(strcmp(ooTable.GFPid,'hb9'));
drd4s = ooTable.vecSum(strcmp(ooTable.GFPid,'drd4'));


% figure, plotSpread({oo_numT, oo_numD, oo_numN, oo_numV, o_numT, o_numD, o_numN, o_numV},'showMM',2)
figure, boxplot([oo_T;oo_D;oo_N;oo_V;o_T;o_D;o_N;o_V;hb9s;drd4s],...
    [ones(length(oo_T),1);ones(length(oo_D),1)*2;ones(length(oo_N),1)*3;ones(length(oo_V),1)*4;...
    ones(length(o_T),1)*5;ones(length(o_D),1)*6;ones(length(o_N),1)*7;ones(length(o_V),1)*8;
    ones(length(hb9s),1)*9; ones(length(drd4s),1)*10],'Notch','on')
axis([0 11 0 1])

figure, violin({oo_N,oo_T,oo_D,oo_V,o_N,o_T,o_D,o_V})
ylim([0 1])
% 
% figure('Name','Horizontal v Vertical VS'), hold
% ecdf([oo_N;oo_T;o_N;o_T])
% ecdf([oo_V;oo_D;o_V;o_D])
% xlim([0 1])

% figure('Name','Hb9 vs Drd4 VS'), hold
% ecdf(hb9s)
% ecdf(drd4s)
% xlim([0 1])

% figure('Name','ON-OFFs VS'), hold
% ecdf(oo_N);
% ecdf(oo_T);
% ecdf(oo_D);
% ecdf(oo_V);
% xlim([0 1])
% 
figure('Name','ONs VS'), hold
ecdf(o_N);
ecdf(o_T);
ecdf(o_D);
ecdf(o_V);
xlim([0 1])

end
