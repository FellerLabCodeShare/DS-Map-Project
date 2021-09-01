function redraw(cellNumb, tuning_fig_pre, tuning_fig_post, dfof_fig_pre, dfof_fig_post , thetas, rhos_ds_pre,rhos_ds_post,...
    ind,rhos1_ds_pre, rhos1_ds_post,rhos2_ds_pre, rhos2_ds_post,rhos3_ds_pre, rhos3_ds_post,vecTheta_ds,vecSum_ds, DSI,...
    DS_roi_dfof_pre, DS_roi_dfof_post, wvf_resp_reordered_pre, wvf_resp_reordered_post, wvf_resp_mean_pre, wvf_resp_mean_post,...
    dir_fig_pre, dir_fig_post,wvf_resp_t1_pre, wvf_resp_t2_pre, wvf_resp_t3_pre, wvf_resp_t1_post, wvf_resp_t2_post, wvf_resp_t3_post,...
    varSum_ds, info_pre, info_post,DSI_ds) 

dfof_llim = -0.2; %Sets the lower limit for dfof in plots below
dfof_ulim = 1.0; %Sets the upper limit for dfof in plots below
shadedColors = [0.8 0.8 0.8];

% Write in some stats for the cell PRE
set(0, 'currentfigure', info_pre);
set(gca,'units','normalized','position',[0,0,1,1])
cla
text(0.1,0.6,'DSI')
text(0.1,0.4,num2str(DSI_ds(cellNumb,1)))
text(0.4,0.6,'VS')
text(0.4,0.4,num2str(vecSum_ds(cellNumb,1)))
text(0.7,0.6,'var')
text(0.7,0.4,num2str(varSum_ds(cellNumb,1)))


% Draw dF/F trace PRE
set(0, 'currentfigure', dfof_fig_pre);
plot(DS_roi_dfof_pre(cellNumb,:))
ylim([dfof_llim dfof_ulim]);
xlim([0,700]);
xlabel('frame')
set(gca,'units','normalized','position',[.035,.025,.96,.95])

% Draw dF/F trace POST
set(0, 'currentfigure', dfof_fig_post);
plot(DS_roi_dfof_post(cellNumb,:))
ylim([dfof_llim dfof_ulim]);
xlim([0,700]);
xlabel('frame')
set(gca,'units','normalized','position',[.035,.025,.96,.95])


% Draw dir trace PRE
set(0, 'currentfigure', dir_fig_pre);
set(dir_fig_pre,'Name',['Cell # ', num2str(ind(cellNumb))],'NumberTitle','off')
cla
patch([13 20 20 13],[-10 -10 10 10],shadedColors, 'LineStyle','none')
hold on
patch([35 42 42 35],[-10 -10 10 10],shadedColors, 'LineStyle','none')
patch([57 64 64 57],[-10 -10 10 10],shadedColors, 'LineStyle','none')
patch([78 85 85 78],[-10 -10 10 10],shadedColors, 'LineStyle','none')
patch([101 108 108 101],[-10 -10 10 10],shadedColors, 'LineStyle','none')
patch([124 131 131 124],[-10 -10 10 10],shadedColors, 'LineStyle','none')
patch([146 153 153 146],[-10 -10 10 10],shadedColors, 'LineStyle','none')
patch([168 175 175 168],[-10 -10 10 10],shadedColors, 'LineStyle','none')
plot(wvf_resp_mean_pre(ind(cellNumb),:),'r', 'LineWidth',2)
plot(wvf_resp_t1_pre(ind(cellNumb),:),'k')
plot(wvf_resp_t2_pre(ind(cellNumb),:),'k')
plot(wvf_resp_t3_pre(ind(cellNumb),:),'k')
hold off

ylim([dfof_llim dfof_ulim]);
xlim([0,176]);
xticks(10:22:176)
xticklabels({'0','45','90','135','180','225','270','315'})
xlabel('angle')
set(gca,'units','normalized','position',[.05,.12,.945,.86])

% Draw dir trace POST
set(0, 'currentfigure', dir_fig_post);
plot(wvf_resp_mean_post(ind(cellNumb),:),'r', 'LineWidth',2)
hold on
plot(wvf_resp_t1_post(ind(cellNumb),:),'k')
plot(wvf_resp_t2_post(ind(cellNumb),:),'k')
plot(wvf_resp_t3_post(ind(cellNumb),:),'k')
hold off
ylim([dfof_llim dfof_ulim]);
xlim([0,176]);
xticks(10:22:176)
xticklabels({'0','45','90','135','180','225','270','315'})
xlabel('angle')
set(gca,'units','normalized','position',[.05,.12,.945,.86])


% Draw PRE-TRAIN tuning curve
set(0, 'currentfigure', tuning_fig_pre);
set(tuning_fig_pre,'Name',['Cell # ', num2str(ind(cellNumb))],'NumberTitle','off')
polarplot(thetas, rhos_ds_pre(cellNumb,:),'k','LineWidth',2)
hold on;
polarplot(thetas, rhos1_ds_pre(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot(thetas, rhos2_ds_pre(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot(thetas, rhos3_ds_pre(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot([0, vecTheta_ds(cellNumb,1)], [0, vecSum_ds(cellNumb,1)]*sum(rhos_ds_pre(cellNumb, 1:8)),'Color', [1,0,0], 'LineWidth',2) %non-normalized vector sum
hold off;
rlim([0 dfof_ulim])

% Draw POST-TRAIN tuning curve
set(0, 'currentfigure', tuning_fig_post);
set(tuning_fig_post,'Name',['Cell # ', num2str(ind(cellNumb))],'NumberTitle','off')
polarplot(thetas, rhos_ds_post(cellNumb,:),'k','LineWidth',2)
hold on;
polarplot(thetas, rhos1_ds_post(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot(thetas, rhos2_ds_post(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot(thetas, rhos3_ds_post(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
polarplot([0, vecTheta_ds(cellNumb,2)], [0, vecSum_ds(cellNumb,2)]*sum(rhos_ds_post(cellNumb, 1:8)),'Color', [1,0,0], 'LineWidth',2) %non-normalized vector sum
hold off;
rlim([0 dfof_ulim])








% Display the tuning curve on right 

% tuning_plot = axes('Units','Normalized','Position',[0.7, 0, 0.3, 0.5], 'XTick',[], 'YTick',[]);
% 
% polarplot(thetas, rhos_ds(cellNumb,:),'k','LineWidth',2)
% hold on;
% polarplot(thetas, rhos1_ds(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
% polarplot(thetas, rhos2_ds(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
% polarplot(thetas, rhos3_ds(cellNumb,:),'Color', [0.2 0.2 0.2],'LineWidth',0.5)
% polarplot([0, vecTheta_ds(cellNumb)], [0, vecSum_ds(cellNumb)]*sum(rhos_ds(cellNumb, 1:8)),'Color', [1,0,0], 'LineWidth',2) %non-normalized vector sum

% tuning_plot.NextPlot = 'replace';



         


end




















