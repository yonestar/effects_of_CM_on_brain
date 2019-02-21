% perm test

% in full (non CV-split) data, find 5% most decreased voxels in Faml, test
% changes in other groups. that is my effect

% permute group labels and repeat 10k times. that gives me p value on my
% effect

% OK: this will ALWAYS show a "buffering" effect, b/c with pure noise
% (e.g., shuffled data), the 5% most decreased voxels in one group will not
% be decreasing in the other groups, hence a buffering effect. 

% So need to do CV within the perm test
% for voxel selection, b/c then with pure noise data, the null should be
% around 0, b/c the observations within one group are independent of the
% others (b/c shuffled). if my effect (g of .25 ish, do repeated CV to get
% a stable estimate) is in the tail of the null, means that ... XXXXXX

clear all

%% load
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
figdir = fullfile(repodir, 'figures');
cd(repodir);

gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));

time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));
delta_relisten = image_math(time2, time1, 'minus');
delta_relisten = apply_mask(delta_relisten, gm_mask);

% divide into the the three groups
load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')
grp = get_var(CM, 'Group');

%% perm test
clear test*

for i=1:100
    i
    grp_shuffle = grp(randperm(length(grp)));

    delta_cm = get_wh_image(delta_relisten, grp_shuffle==1);
    delta_oxy = get_wh_image(delta_relisten, grp_shuffle==2);
    delta_faml = get_wh_image(delta_relisten, grp_shuffle==3);

    
    %% find 5% most decreased voxels in "Faml" set (in quotes b/c shuffled)
    t = ttest(delta_faml);
    decreases5prc = threshold(t, [prctile(t.dat, 5) Inf], 'raw-outside');
    decreases5prc.dat = decreases5prc.sig; % make into a "mask"

    %decreases5prc = threshold(t, .01, 'unc');
    %figure; montage(decreases5prc)
    
    %% extract ROI vals from other "groups" (in quotes b/c shuffled)
    
    tmp = extract_roi_averages(delta_relisten, decreases5prc);
    test_cm(:,i) = tmp.dat(grp_shuffle==1);
    test_oxy(:,i) = tmp.dat(grp_shuffle==2);
    test_faml(:,i) = tmp.dat(grp_shuffle==3);
    
    
end % end perm iteration

save(fullfile(repodir, 'results', 'permtestresults.mat'), 'test_cm','test_oxy','test_faml');

%% look at null distribution 

figure; hist(mean(test_cm) - mean(test_faml)), hold on


%% compute and plot the real effect

delta_cm = get_wh_image(delta_relisten, grp==1);
delta_oxy = get_wh_image(delta_relisten, grp==2);
delta_faml = get_wh_image(delta_relisten, grp==3);

% find 5% most decreased voxels in Faml
t = ttest(delta_faml);
decreases5prc = threshold(t, [prctile(t.dat, 5) Inf], 'raw-outside');
decreases5prc.dat = decreases5prc.sig; % make into a "mask"

% extract ROI vals from other groups
tmp = extract_roi_averages(delta_relisten, decreases5prc);
cm_vs_faml_buff = mean(tmp.dat(grp==1)) - mean(tmp.dat(grp==3))
oxy_vs_faml_buff = mean(tmp.dat(grp==2)) - mean(tmp.dat(grp==3))

vline(cm_vs_faml_buff)


%% save output

%save(fullfile(repodir,'data','cv_buffering_5prc_splithalf.mat'), 'cm_vs_faml', 'oxy_vs_faml', 'cm_vs_oxy')

% plot distributions

figure
violinplot({test_cm test_oxy test_faml},10)

% test

g = mes(test_cm, test_faml, 'hedgesg', 'nBoot', 1000)
%g2 = mes(oxy_vs_faml(:), 0, 'hedgesg', 'nBoot', 1000)

[~,p]=ttest2(test_cm, test_faml)