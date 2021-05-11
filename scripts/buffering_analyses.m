% in CV loop, find 5% most decreased voxels in Faml, test
% changes in held out data from all 3 groups. that is my true effect.
% permute group labels and repeat 10k times. that gives me p value on my
% effect


%% set up
clear all

repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
figdir = fullfile(repodir, 'figures');
cd(repodir); cd scripts

gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));

time1 = fmri_data(filenames(fullfile(repodir,'data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile(repodir,'data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));
delta_relisten = image_math(time2, time1, 'minus');
delta_relisten = apply_mask(delta_relisten, gm_mask);

load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')
grp = get_var(CM, 'Group');


%% do repeated CV to get a point estimate for the true effect
for i=1:10
    i
    [test_cm(i,:), test_oxy(i,:), test_faml(i,:), ~] = cv_buffering_test(delta_relisten, grp, 5);
    g = mes(test_cm(i,:), test_faml(i,:), 'hedgesg'); 
    g_cm_vs_faml(i) = g.hedgesg;
end

%% get a null distribution: permutation test with shuffled group labels
for i=1:10
    i
    grp_shuffle = grp(randperm(length(grp)));
    [test_cm(i,:), test_oxy(i,:), test_faml(i,:), decreases{i}] = cv_buffering_test(delta_relisten, grp_shuffle, 5);
    g = mes(test_cm(i,:), test_faml(i,:), 'hedgesg'); 
    g_cm_vs_faml_null(i) = g.hedgesg;
end

%% save output

save(fullfile(repodir,'data','cv_buffering_perm_test.mat'), 'g_cm_vs_faml', 'g_cm_vs_faml_null')

%% visualize the distribution of repeated CV effects, just for curiousity and checking
figure; histogram(g_cm_vs_faml)


%% plot the null distribution
figure; histogram(g_cm_vs_faml_null)
vline(mean(g_cm_vs_faml))