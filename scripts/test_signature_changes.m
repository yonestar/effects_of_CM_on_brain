%% set up and load data
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
cd(repodir);

% load results/group_level/CM_vs_Faml_robust.mat
% dat = get_wh_image(out.t,1);

% load subject level data
time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

% time2 - time1
delta_relisten = image_math(time2, time1, 'minus'); 

% behav data
load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat')) % load behav data
grp = get_var(CM, 'Group');
%% load EC ED models

path = '/Users/yoni/Repositories/Neuroimaging_Pattern_Masks/Multivariate_signature_patterns/2017_Ashar_care_distress';
ec = fmri_data(fullfile(path, 'Ashar_2017_empathic_care_marker.nii'));
ed = fmri_data(fullfile(path, 'Ashar_2017_empathic_distress_marker.nii'));


%% apply patterns: EC

delta_ec = apply_mask(delta_relisten, ec, 'pattern_expression');
figure; boxchart(grp, delta_ec)
[~, p, ~, stats] = ttest2(delta_ec(grp==1), delta_ec(grp==3))

%% apply patterns: ED

delta_ed = apply_mask(delta_relisten, ed, 'pattern_expression');
figure; boxchart(grp, delta_ed)
[~, p, ~, stats] = ttest2(delta_ed(grp==1), delta_ed(grp==3))



%%  OLD
% before running this script, must run prep_4_apply_signatures_and_save.m
% this saves signature response values in the DAT object.

basedir = '/Users/yoni/Repositories/effects_of_CM_on_brain';

load(fullfile(basedir, 'results', 'image_names_and_setup.mat'));

distress = DAT.SIG_contrasts.raw.cosine_sim.Empathic_Dist;
care = DAT.SIG_contrasts.raw.cosine_sim.Empathic_Care;

group = DAT.BETWEENPERSON.group;

d{1} = distress.pre_to_post_change