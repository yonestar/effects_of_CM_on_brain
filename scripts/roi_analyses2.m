
% roi anaylses

projdir = '/dreamio3/wagerlab/labdata/projects/LKM_fmri_longitudinal/'; 
datadir = '/dreamio3/wagerlab/labdata/data/LKM/'; 

behavdatfile = fullfile(projdir,'CM_data_2015-07-29_16_16.mat');
load(behavdatfile)

projdir2 = fullfile(projdir, 'Analyses');

warning off

%%
wh = CM.wh_keep.t2; 
group = get_var(CM, 'Group', wh);
deltadon = get_var(CM, 'deltadon', wh);
deltaFAS = get_var(CM, 'deltaFAS', wh);

endind = size(CM.Subj_Level.data,2);

%% legend
fprintf('Red is CM, Blue is Oxy, Black is Faml\n');

%% load anamtomically based ROI masks from folder
roidir = fullfile(projdir2, 'roi_analysis', 'roimasks');%, 'ns_masks');
roi_masks = dir(fullfile(roidir, '*nii')); 

% only use one mask
%roi_masks = fullfile(projdir, 'roi_analysis', 'roimasks','NAcc_5perc_and_up.nii');

% data from contrast of interest
subjdat_t1 = fmri_data(canlab_list_files(projdir2, 'first_level', 'time1and2_no_pmods', CM.Subj_Level.id(wh), 'con_0003.img'));
subjdat_t2 = fmri_data(canlab_list_files(projdir2, 'first_level', 'time1and2_no_pmods', CM.Subj_Level.id(wh), 'con_0007.img'));

%t1dat = fmri_data(canlab_list_files(projdir, 'first_level', 'tenddistdon_idiographic_timecourses', CM.Subj_Level.id(wh), 'con_0004.img')); %don by don amounts
i=1; clear size;
%% get ventricle activity for this contrast
% ventricle = image_vector('image_names', fullfile(projdir, 'anat_masks', 'ventricles.img'));
% ventricle.dat = ventricle.dat > 0;
% v1 = extract_roi_averages(subjdat_t1, ventricle);
% v2 = extract_roi_averages(subjdat_t2, ventricle);

%% for each mask
for i=1:length(roi_masks)

    mask = image_vector('image_names', fullfile(roidir, roi_masks(i).name));
    %mask = image_vector('image_names', fullfile(roi_masks));
    
    %binarize and only take positive
    mask.dat = mask.dat > 0;
    orthviews(mask, 'largest_region');
    snapnow
    
    %% extract time 1 and 2 data from ROIs
    cl1 = extract_roi_averages(subjdat_t1, mask); 
    cl2 = extract_roi_averages(subjdat_t2, mask);
    
    % ROI change minus ventricle change
    %cl1.dat = cl1.dat - v1.dat;  cl2.dat = cl2.dat - v2.dat;
    clear roiact;
    roiact(:,:,i) = [cl1.dat cl2.dat (cl2.dat - cl1.dat)]; 

    %add to CM temporarily
    CM.Subj_Level.data(:,endind+1:endind+3) = zeros(size(CM.Subj_Level.data,1),3);
    CM.Subj_Level.data(wh,endind+1:endind+3) = roiact(:,:,i);
    CM.Subj_Level.names(endind+1:endind+3) = {'roi_t1', 'roi_t2', 'roi_delta'};
    
    name = roi_masks(i).name;
    %[~,name]=fileparts(roi_masks);  % with only one mask
     
    % plot
    h=create_figure(name, 2, 2);     title(name);
    set(h, 'Position', [0 0 1200 700]);

    
    % PANE 1: line plot
    subplot(2,2,1);
    lineplot_columns(roiact(group==1,1:2,i), 'color', 'r');
    lineplot_columns(roiact(group==2,1:2,i), 'color', 'b');
    lineplot_columns(roiact(group==3,1:2,i), 'color', 'k');
    xlim([.8 2.2]);
    
    % PANE 2: bar graph
    subplot(2,2,2);
    plot_var(CM, 'roi_delta', 'subjtype', 'Group', 'wh_keep', wh, '95CI', 'nofig');
    
    % PANES 3, 4:  scatter plots with deltaFAS and deltadon
    subplot(2,2,3); xlabel('delta brain'); ylabel('delta FAS');
    scatterplot(CM, 'roi_delta', 'deltaFAS', 'subjtype', 'Group', 'wh_keep', wh, 'nofig');
   
    subplot(2,2,4); xlabel('delta brain'); ylabel('delta donation');
    scatterplot(CM, 'roi_delta', 'deltadon', 'subjtype', 'Group', 'wh_keep', wh, 'nofig');
        
    %subplot(2,3,4);
    %scatterplot(CM, 'roi_delta', 'wd1', 'subjtype', 'Group', 'wh_keep', wh, 'nofig');
    %scatterplot(CM, 'roi_delta', 'wd1', 'wh_keep', wh, 'nofig', 'dorobust');
    %plot_correlation_samefig(roiact(:,3), get_var(CM,'deltadon',wh),[],[],0,1)
    %xlabel('delta brain'); ylabel('delta donation');

    % tenderness
%     subplot(2,3,5);
%     scatterplot(CM, 'roi_delta', 'DeltaTend_mean', 'subjtype', 'Group', 'wh_keep', wh, 'nofig');
%     xlabel('delta brain'); ylabel('delta tend');
%     xlim([-.3 .6])
%     
%     subplot(2,3,6);
%     plot_correlation_samefig(roiact(:,3), get_var(CM,'DeltaTend_mean',wh),[],[],0,1)
%     xlabel('delta brain'); ylabel('deltatenderness'); xlim([-.8 1])
%     
%     
    
    %% capture output
    fprintf('\n\n\n\nName of mask:  %s.\n  Note that only positive activations were taken for the TICS maps!', roi_masks(i).name);
    snapnow
    close all
     
    %t ttests
    [~,p, ci, stats] = ttest2(CM, 'roi_delta', wh & CM.wh_keep.med, wh & ~CM.wh_keep.med, 'noverbose');
    fprintf('\n\n%s: CM vs. other, t(%CM) = %3.2f\n', roi_masks(i).name, stats.df, p);
    [~,p, ci, stats] = ttest2(CM, 'roi_delta', wh & CM.wh_keep.med, wh & CM.wh_keep.faml, 'noverbose');
    fprintf('\n\n%s: CM vs. faml, t(%CM) = %3.2f\n\n\n', roi_masks(i).name, stats.df, p);

end

%{

%% just do scatter plots for a given ROI
roidir = fullfile(projdir, 'roi_analysis', 'roimasks');
mask = image_vector('image_names', fullfile(roidir, 'rsgACC_p_01.nii'));
orthviews(mask, 'largest_region')
%%
subjdat_delta = fmri_data(canlab_list_files(projdir, 'first_level', 'time1and2_no_pmods', CM.Subj_Level.id, 'con_0014.img'));
cl = extract_roi_averages(subjdat_delta, mask, 'contiguous_regions');
%% add to CM temporarily
endind = size(CM.Subj_Level.data,2);
CM.Subj_Level.data(:,endind+1) = zeros(size(CM.Subj_Level.data,1),1);
CM.Subj_Level.data(:,endind+1) = cl.dat;
CM.Subj_Level.names(endind+1) = {'roi_delta'};
%% scatter plots
create_figure('ROI analysis, delta data only');
scatterplot(CM, 'roi_delta', 'ccd1', 'subjtype', 'Group')%, 'wh_keep', CM.wh_keep.allmed);
%scatterplot(CM, 'roi_delta', 'ccd1', 'wh_keep', wh, 'nofig', 'dorobust');
xlabel('delta brain'); ylabel('delta behav');
create_figure('robust')
wh=CM.wh_keep.all;
plot_correlation_samefig(cl.dat(wh), get_var(CM, 'ccd1',wh), [],[],0,1)



%% are Ss high in one ROI high in others too? hard to visually see, but mean corr significantly above zero

roi_means = squeeze(roiact(:,3,:))
create_figure('rois delta by subj'); plot(roi_means(1:10,:)')%lineplot_columns(roi_means)
mean(corr(roi_means)) % mean corr is clearly above zero.
ttest(corr(roi_means))

%% look at change in the ventricles
% ideally, get an avg gm and wm mask for my Ss (use
% canlab_create_wm_ventricle) and remove that from
% canlab_canonical_ventricle.  
anats = canlab_list_files(datadir, 'Imaging', CM.Subj_Level.id(wh), 'Structural', 'SPGR', 'wmprage_avg.nii');
gms = canlab_list_files(datadir, 'Imaging', CM.Subj_Level.id(wh), 'Structural', 'SPGR', 'wc1mprage_avg.nii');
wms = canlab_list_files(datadir, 'Imaging', CM.Subj_Level.id(wh), 'Structural', 'SPGR', 'wc2mprage_avg.nii');

canlab_combine_structurals_luka(fullfile(projdir, 'anat_masks', 'avg_wc1_gm.nii'), gms);
canlab_combine_structurals_luka(fullfile(projdir, 'anat_masks', 'avg_wc2_wm.nii'), wms);
canlab_combine_structurals_luka(fullfile(projdir, 'anat_masks', 'avg_wmprage.nii'), anats);

%%
canlab_create_wm_ventricle_masks( fullfile(projdir, 'anat_masks', 'avg_wc2_wm.nii') , fullfile(projdir, 'anat_masks','avg_wc1_gm.nii'), .3 )



%% look at ROIs from main effect of time in donation epoch

%get ROIs from main effect of time on don epoch
img = image_vector('image_names', fullfile(basedir, 'Analyses', 'first_level', 'no_timecourses_time1and2', 'deltaFAS', 'robust0001', 'rob_p_0001.img'))
%img.p = img.dat;
mask_rs = fullfile(basedir, 'Analyses', 'scalped_avg152T1_graymatter_smoothed_resampleCM.img');

img = apply_mask(img, mask_rs);
img = threshold(img, [0 .05], 'raw-between', 'k', 900);
orthviews(img)

don_mask = img; %cl 1 is large posterior region, 2 = lVMPFC, 3 = right temporal and STS area, 4 = left STS, 5 =  
don_mask.dat = img.dat > 0; 

  %}  