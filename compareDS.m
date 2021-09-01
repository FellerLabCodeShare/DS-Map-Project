%% Compare DS pre and post Training

%Loading the calcium movies
disp('pick the pre-training calcium imaging file'); 
[movie_preTrain, path_name] = uigetfile([path, 'experiment', '\*.tif']);

disp('pick the post-training calcium imaging file'); 
[movie_postTrain, path_name] = uigetfile([path, 'experiment', '\*.tif']);

%Loading the mask
disp('pick the calcium ROI mask file'); 
[roi_mask_file, path_name] = uigetfile([path, 'experiment', '\*.tif']);

%Loading the text files
disp('pick the pre-training text file with DS directions')
[textF_preTrain, path_name_text] = uigetfile([path, 'experiment', '\*.txt']);

disp('pick the post-training text file with DS directions')
[textF_postTrain, path_name_text] = uigetfile([path, 'experiment', '\*.txt']);

textFile_preTrain = readtable(textF_preTrain);
textFile_postTrain = readtable(textF_postTrain);

textFileArray_preTrain = table2array(textFile_preTrain);
textFileArray_postTrain = table2array(textFile_postTrain);


%Loading the movies
n_frames = numel(imfinfo(movie_preTrain));
[height, width] = size(imread(movie_preTrain,1));
movie_pre = zeros(height, width, n_frames);
movie_post = zeros(height, width, n_frames);
for k = 1:n_frames
    movie_pre(:,:,k) = imread(movie_preTrain, k);
    movie_post(:,:,k) = imread(movie_postTrain, k);
end
%End loading the movies

%Loading the mask
roi_mask = imread(roi_mask_file);
countmask=bwlabel(roi_mask);
%End loading the mask

%% Calibrate the directions on the text file
r0 = 0;
r45 = 45;
r90 = 90;
r135 = 135;
r180 = 180;
r225 = 225;
r270 = 270;
r315 = 315;





% r0 = 18;
% r45 = 63;
% r90 = 108;
% r135 = 153;
% r180 = 198;
% r225 = 243;
% r270 = 288;
% r315 = 333;

temp = textFileArray_preTrain(:,1);
if min(temp) > 0
    temp(temp == r0) = 0;
    temp(temp == r45) = 45;
    temp(temp == r90) = 90;
    temp(temp == r135) = 135;
    temp(temp == r180) = 180;
    temp(temp == r225) = 225;
    temp(temp == r270) = 270;
    temp(temp == r315) = 315;
    textFileArray_preTrain(:,1) = temp;
    textFileArray_postTrain(:,1) = temp;
end


%% use ROI mask to find z profile at each cell

for i=1:n_frames
    partprops_pre=regionprops(countmask,movie_pre(:,:,i),'Area','MeanIntensity'); 
    partprops_post=regionprops(countmask,movie_post(:,:,i),'Area','MeanIntensity');
    for j = 1:max(max(countmask))
        roiInt_pre(j,i) = partprops_pre(j).MeanIntensity;
        roiInt_post(j,i) = partprops_post(j).MeanIntensity;
    end
end

[num_cells, num_frames] = size(roiInt_pre);


% load('backF');  %This set of code is to remove bleedthrough
% backF_exp = repmat(backgroundF',[num_cells,1]);
% roiInt2 = roiInt-backF_exp;
% roiInt = roiInt2;


%% Figures for looking at cell responses across all stims
% figure, imagesc(roiInt_pre) % plot the whole pre movie as an image of dFoF
% figure, imagesc(roiInt_post) % plot the whole post movie as an image of dFoF

%% Find max intensity at every bar presentation (not discriminating between ON and OFF)

frameRate = 1.48;  %1.48

%Pre
barResp_pre = zeros(num_cells,24);
wvf_resp_pre = zeros(num_cells, 22,24); % waveforms of all responses

%Post
barResp_post = zeros(num_cells,24);
wvf_resp_post = zeros(num_cells, 22,24); % waveforms of all responses

for i = 1:24
    %Pre
    stimFrameStart_pre = floor(textFileArray_preTrain(i,7)*frameRate);
    stimFrameEnd_pre = ceil(textFileArray_preTrain(i,8)*frameRate);
    barResp_pre(:,i) = max(roiInt_pre(:,stimFrameStart_pre:stimFrameEnd_pre),[],2);
    
    wvf_resp_pre(:,:,i) = roiInt_pre(:,stimFrameStart_pre-3:stimFrameStart_pre+18);
    
    %Post
    stimFrameStart_post = floor(textFileArray_postTrain(i,7)*frameRate);
    stimFrameEnd_post = ceil(textFileArray_postTrain(i,8)*frameRate);
    barResp_post(:,i) = max(roiInt_post(:,stimFrameStart_post:stimFrameEnd_post),[],2);
    
    wvf_resp_post(:,:,i) = roiInt_post(:,stimFrameStart_post-3:stimFrameStart_post+18);    
end

%Pre
ind0_pre = find(textFileArray_preTrain(:,1) == 0);
ind45_pre = find(textFileArray_preTrain(:,1) == 45);
ind90_pre = find(textFileArray_preTrain(:,1) == 90);
ind135_pre = find(textFileArray_preTrain(:,1) == 135);
ind180_pre = find(textFileArray_preTrain(:,1) == 180);
ind225_pre = find(textFileArray_preTrain(:,1) == 225);
ind270_pre = find(textFileArray_preTrain(:,1) == 270);
ind315_pre = find(textFileArray_preTrain(:,1) == 315);

ind0_post = find(textFileArray_postTrain(:,1) == 0);
ind45_post = find(textFileArray_postTrain(:,1) == 45);
ind90_post = find(textFileArray_postTrain(:,1) == 90);
ind135_post = find(textFileArray_postTrain(:,1) == 135);
ind180_post = find(textFileArray_postTrain(:,1) == 180);
ind225_post = find(textFileArray_postTrain(:,1) == 225);
ind270_post = find(textFileArray_postTrain(:,1) == 270);
ind315_post = find(textFileArray_postTrain(:,1) == 315);


%% now tuning curves

%Zeros start
DSI = zeros(num_cells,2);
vecSum = zeros(num_cells,2);
vecTheta = zeros(num_cells,2);
maxDF = zeros(num_cells,2);
varSum = zeros(num_cells,2);

rhos_all_pre = zeros(num_cells, 9);
rhos1_all_pre = zeros(num_cells, 9);
rhos2_all_pre = zeros(num_cells, 9);
rhos3_all_pre = zeros(num_cells, 9);
rhos_var_pre = zeros(num_cells, 9);

rhos_all_post = zeros(num_cells, 9);
rhos1_all_post = zeros(num_cells, 9);
rhos2_all_post = zeros(num_cells, 9);
rhos3_all_post = zeros(num_cells, 9);
rhos_var_post = zeros(num_cells, 9);
%Zeroes end

% temp_barResp = barResp_pre;
% barResp_pre = barResp_pre - 1; %Uncomment this ifyou want to remove back
% barResp_pre(barResp_pre<0) = 0;


thetas = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 0]; %Theta values for 0 45 90 135 180 225 270 315


for i = 1:num_cells
    
    rhos_pre = [mean(barResp_pre(i,ind0_pre)), mean(barResp_pre(i,ind45_pre)), mean(barResp_pre(i,ind90_pre)), mean(barResp_pre(i,ind135_pre)), mean(barResp_pre(i,ind180_pre)), mean(barResp_pre(i,ind225_pre)), mean(barResp_pre(i,ind270_pre)), mean(barResp_pre(i,ind315_pre)), mean(barResp_pre(i,ind0_pre))];
    rhos1_pre = [barResp_pre(i,ind0_pre(1)), barResp_pre(i,ind45_pre(1)), barResp_pre(i,ind90_pre(1)), barResp_pre(i,ind135_pre(1)), barResp_pre(i,ind180_pre(1)), barResp_pre(i,ind225_pre(1)), barResp_pre(i,ind270_pre(1)), barResp_pre(i,ind315_pre(1)), barResp_pre(i,ind0_pre(1))];
    rhos2_pre = [barResp_pre(i,ind0_pre(2)), barResp_pre(i,ind45_pre(2)), barResp_pre(i,ind90_pre(2)), barResp_pre(i,ind135_pre(2)), barResp_pre(i,ind180_pre(2)), barResp_pre(i,ind225_pre(2)), barResp_pre(i,ind270_pre(2)), barResp_pre(i,ind315_pre(2)), barResp_pre(i,ind0_pre(2))];
    rhos3_pre = [barResp_pre(i,ind0_pre(3)), barResp_pre(i,ind45_pre(3)), barResp_pre(i,ind90_pre(3)), barResp_pre(i,ind135_pre(3)), barResp_pre(i,ind180_pre(3)), barResp_pre(i,ind225_pre(3)), barResp_pre(i,ind270_pre(3)), barResp_pre(i,ind315_pre(3)), barResp_pre(i,ind0_pre(3))];
    rhos_var_pre(i,:) = [var(barResp_pre(i,ind0_pre)), var(barResp_pre(i,ind45_pre)), var(barResp_pre(i,ind90_pre)), var(barResp_pre(i,ind135_pre)), var(barResp_pre(i,ind180_pre)), var(barResp_pre(i,ind225_pre)), var(barResp_pre(i,ind270_pre)), var(barResp_pre(i,ind315_pre)), var(barResp_pre(i,ind0_pre))];
    varSum(i,1) = sum(rhos_var_pre(i,:));
        
    rhos_post = [mean(barResp_post(i,ind0_post)), mean(barResp_post(i,ind45_post)), mean(barResp_post(i,ind90_post)), mean(barResp_post(i,ind135_post)), mean(barResp_post(i,ind180_post)), mean(barResp_post(i,ind225_post)), mean(barResp_post(i,ind270_post)), mean(barResp_post(i,ind315_post)), mean(barResp_post(i,ind0_post))];
    rhos1_post = [barResp_post(i,ind0_post(1)), barResp_post(i,ind45_post(1)), barResp_post(i,ind90_post(1)), barResp_post(i,ind135_post(1)), barResp_post(i,ind180_post(1)), barResp_post(i,ind225_post(1)), barResp_post(i,ind270_post(1)), barResp_post(i,ind315_post(1)), barResp_post(i,ind0_post(1))];
    rhos2_post = [barResp_post(i,ind0_post(2)), barResp_post(i,ind45_post(2)), barResp_post(i,ind90_post(2)), barResp_post(i,ind135_post(2)), barResp_post(i,ind180_post(2)), barResp_post(i,ind225_post(2)), barResp_post(i,ind270_post(2)), barResp_post(i,ind315_post(2)), barResp_post(i,ind0_post(2))];
    rhos3_post = [barResp_post(i,ind0_post(3)), barResp_post(i,ind45_post(3)), barResp_post(i,ind90_post(3)), barResp_post(i,ind135_post(3)), barResp_post(i,ind180_post(3)), barResp_post(i,ind225_post(3)), barResp_post(i,ind270_post(3)), barResp_post(i,ind315_post(3)), barResp_post(i,ind0_post(3))];
    rhos_var_post(i,:) = [var(barResp_post(i,ind0_post)), var(barResp_post(i,ind45_post)), var(barResp_post(i,ind90_post)), var(barResp_post(i,ind135_post)), var(barResp_post(i,ind180_post)), var(barResp_post(i,ind225_post)), var(barResp_post(i,ind270_post)), var(barResp_post(i,ind315_post)), var(barResp_post(i,ind0_post))];
    varSum(i,2) = sum(rhos_var_post(i,:));
  
    rhosNorm_pre = rhos_pre/sum(rhos_pre(1:8));
    [x_pre,y_pre] = pol2cart(thetas,rhosNorm_pre);
    vectSumX_pre = sum(x_pre(1:8));
    vectSumY_pre = sum(y_pre(1:8));
    [the_pre, rho_pre] = cart2pol(vectSumX_pre,vectSumY_pre);
    if the_pre <0
        the_pre = 2*pi+the_pre;
    end
    
    rhosNorm_post = rhos_post/sum(rhos_post(1:8));
    [x_post,y_post] = pol2cart(thetas,rhosNorm_post);
    vectSumX_post = sum(x_post(1:8));
    vectSumY_post = sum(y_post(1:8));
    [the_post, rho_post] = cart2pol(vectSumX_post,vectSumY_post);
    if the_post <0
        the_post = 2*pi+the_post;
    end
    
    %what is the pref dir?
    
    
    if the_pre > 5.8905
        prefInd_pre = 1;
    else
        [temp,prefInd_pre] = min(abs(thetas-the_pre));
    end
    prefDir = thetas(prefInd_pre);
    nullInd_pre = prefInd_pre - 4;
    if nullInd_pre < 1
        nullInd_pre = 8+nullInd_pre;
    end
    
    DSI_temp_pre = (rhos_pre(prefInd_pre)-rhos_pre(nullInd_pre))/(rhos_pre(prefInd_pre)+rhos_pre(nullInd_pre));
    
    if the_post > 5.8905
        prefInd_post = 1;
    else
        [temp,prefInd_post] = min(abs(thetas-the_post));
    end
    prefDir = thetas(prefInd_post);
    nullInd_post = prefInd_post - 4;
    if nullInd_post < 1
        nullInd_post = 8+nullInd_post;
    end
    
    DSI_temp_post = (rhos_post(prefInd_post)-rhos_post(nullInd_post))/(rhos_post(prefInd_post)+rhos_post(nullInd_post));
    
    
    %Store for cells
    DSI(i,1) = DSI_temp_pre;
    DSI(i,2) = DSI_temp_post;
    vecSum(i,1) = rho_pre;
    vecSum(i,2) = rho_post;
    vecTheta(i,1) = the_pre;
    vecTheta(i,2) = the_post;
    maxDF(i,1) = max(rhos_pre); %What was the max DF
    maxDF(i,2) = max(rhos_post); %What was the max DF
    
    rhos_all_pre(i,:) = rhos_pre;
    rhos1_all_pre(i,:) = rhos1_pre;
    rhos2_all_pre(i,:) = rhos2_pre;
    rhos3_all_pre(i,:) = rhos3_pre;  
    
    rhos_all_post(i,:) = rhos_post;
    rhos1_all_post(i,:) = rhos1_post;
    rhos2_all_post(i,:) = rhos2_post;
    rhos3_all_post(i,:) = rhos3_post;
    %End store for cells
    

end


figure
subplot(1,3,1)
plot(vecSum(:,1), DSI(:,1),'k.');
hold
plot(vecSum(:,2), DSI(:,2),'r.');
xlabel('Norm Vector Sum')
ylabel('DSI')
subplot(1,3,2)
plot(DSI(:,1),DSI(:,2),'k.')
axis([-.2 0.8 -.2 0.8])
subplot(1,3,3)
plot(vecSum(:,1),vecSum(:,2),'k.')



%% Interactive plot

%Include only cells with following thresholds
DSI_thresh = -1;
vecSum_thresh = 0;
maxDF_thresh = 0;

ind = find(DSI(:,1)>DSI_thresh & vecSum(:,1)>vecSum_thresh & maxDF(:,1) > maxDF_thresh);


% Use wvf_resp to make plots of mean waveforms for each dir for each cell

reordered_ind_pre = [ind0_pre;ind45_pre;ind90_pre;ind135_pre;ind180_pre;ind225_pre;ind270_pre;ind315_pre];
reordered_ind_post = [ind0_post;ind45_post;ind90_post;ind135_post;ind180_post;ind225_post;ind270_post;ind315_post];

wvf_resp_reordered_pre = wvf_resp_pre(:,:,reordered_ind_pre);
wvf_resp_reordered_post = wvf_resp_post(:,:,reordered_ind_post);

wvf_resp_mean_pre = [mean(wvf_resp_reordered_pre(:,:,1:3),3),mean(wvf_resp_reordered_pre(:,:,4:6),3),mean(wvf_resp_reordered_pre(:,:,7:9),3),mean(wvf_resp_reordered_pre(:,:,10:12),3),mean(wvf_resp_reordered_pre(:,:,13:15),3),mean(wvf_resp_reordered_pre(:,:,16:18),3),mean(wvf_resp_reordered_pre(:,:,19:21),3),mean(wvf_resp_reordered_pre(:,:,22:24),3)];
wvf_resp_mean_post = [mean(wvf_resp_reordered_post(:,:,1:3),3),mean(wvf_resp_reordered_post(:,:,4:6),3),mean(wvf_resp_reordered_post(:,:,7:9),3),mean(wvf_resp_reordered_post(:,:,10:12),3),mean(wvf_resp_reordered_post(:,:,13:15),3),mean(wvf_resp_reordered_post(:,:,16:18),3),mean(wvf_resp_reordered_post(:,:,19:21),3),mean(wvf_resp_reordered_post(:,:,22:24),3)];

wvf_resp_t1_pre = []; %These 6 lines just creates the variables
wvf_resp_t2_pre = [];
wvf_resp_t3_pre = [];
wvf_resp_t1_post = [];
wvf_resp_t2_post = [];
wvf_resp_t3_post = [];
for i = 1:8 %This for loop will create wvf for each trial (3 trials here)
    wvf_resp_t1_pre = [wvf_resp_t1_pre, wvf_resp_reordered_pre(:,:,i*3-2)];
    wvf_resp_t2_pre = [wvf_resp_t2_pre, wvf_resp_reordered_pre(:,:,i*3-1)];
    wvf_resp_t3_pre = [wvf_resp_t3_pre, wvf_resp_reordered_pre(:,:,i*3)];
    
    wvf_resp_t1_post = [wvf_resp_t1_post, wvf_resp_reordered_post(:,:,i*3-2)];
    wvf_resp_t2_post = [wvf_resp_t2_post, wvf_resp_reordered_post(:,:,i*3-1)];
    wvf_resp_t3_post = [wvf_resp_t3_post, wvf_resp_reordered_post(:,:,i*3)];
end


% Figures to scroll through cells

DS_roi_dfof_pre = roiInt_pre(ind,:);
rhos_ds_pre = rhos_all_pre(ind, :);
rhos1_ds_pre = rhos1_all_pre(ind, :);
rhos2_ds_pre = rhos2_all_pre(ind, :);
rhos3_ds_pre = rhos3_all_pre(ind, :);

DS_roi_dfof_post = roiInt_post(ind,:);
rhos_ds_post = rhos_all_post(ind, :);
rhos1_ds_post = rhos1_all_post(ind, :);
rhos2_ds_post = rhos2_all_post(ind, :);
rhos3_ds_post = rhos3_all_post(ind, :);

DSI_ds = DSI(ind,:);
vecSum_ds = vecSum(ind,:);
vecTheta_ds = vecTheta(ind,:);
varSum_ds = varSum(ind,:);

cell_numb = 1;


% Create figure windows and save the handles 
info_pre = figure; 
set(gcf, 'Position', [10   100   601   100]);
dfof_fig_pre = figure; 
set(gcf, 'Position', [10   200   601   100]);
dir_fig_pre = figure; 
set(gcf, 'Position', [10   380   400   200]);
tuning_fig_pre = figure; 
set(gcf, 'Position', [411   380   200   200]);

info_post = figure; 
set(gcf, 'Position', [620   100   601   100]);
dfof_fig_post = figure; 
set(gcf, 'Position', [620   200   601   100]);
dir_fig_post = figure; 
set(gcf, 'Position', [620   380   400   200]);
tuning_fig_post = figure; 
set(gcf, 'Position', [1021   380   200   200]);


% Create wrapper function for the redraw function called by video_fig so
% that I can pass in 'constant' arguments in addition to updating frame
% argument
redraw_wrapper = @(cellNumb) redraw(cellNumb, tuning_fig_pre, tuning_fig_post, dfof_fig_pre, dfof_fig_post, thetas, rhos_ds_pre,rhos_ds_post,...
    ind,rhos1_ds_pre, rhos1_ds_post,rhos2_ds_pre, rhos2_ds_post,rhos3_ds_pre, rhos3_ds_post,vecTheta_ds,vecSum_ds, DSI, DS_roi_dfof_pre,...
    DS_roi_dfof_post, wvf_resp_reordered_pre, wvf_resp_reordered_post, wvf_resp_mean_pre, wvf_resp_mean_post, dir_fig_pre, dir_fig_post,...
    wvf_resp_t1_pre, wvf_resp_t2_pre, wvf_resp_t3_pre, wvf_resp_t1_post, wvf_resp_t2_post, wvf_resp_t3_post, varSum_ds, info_pre, info_post,...
    DSI_ds); 

% Initialize to the start of the movie
redraw(1, tuning_fig_pre, tuning_fig_post, dfof_fig_pre, dfof_fig_post, thetas, rhos_ds_pre, rhos_ds_post,ind,rhos1_ds_pre, rhos1_ds_post,...
    rhos2_ds_pre, rhos2_ds_post,rhos3_ds_pre, rhos3_ds_post,vecTheta_ds,vecSum_ds, DSI, DS_roi_dfof_pre, DS_roi_dfof_post,...
    wvf_resp_reordered_pre, wvf_resp_reordered_post, wvf_resp_mean_pre, wvf_resp_mean_post, dir_fig_pre, dir_fig_post,...
    wvf_resp_t1_pre, wvf_resp_t2_pre, wvf_resp_t3_pre, wvf_resp_t1_post, wvf_resp_t2_post, wvf_resp_t3_post, varSum_ds, info_pre, info_post,...
    DSI_ds);

% Hand control over to the user
[scroller, axes_handle, scroll_bar_handles, scroll_func] = videofig(length(ind), redraw_wrapper);