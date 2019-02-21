clear all

%% load data and other set up
 
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
outdir = fullfile(repodir, 'results', 'group_level');
cd(repodir);

time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

% time2 - time1
delta_relisten = image_math(time2, time1, 'minus'); 

% load group assignments and behavioral data from trial
% this data is sorted by subject num in alphanumeric order
% so it maps on to the imaging data one-to-one, in order
load(fullfile('data', 'CM_data_2015-07-29_16_16.mat'))  % canlab dataset object.  loads 'CM' var

% drop the 13 Ps who were not randomized (see manuscript). they don't have
% brain data here, so don't need to drop brain data.
CM.Subj_Level.data(1:13, :) = [];
CM.Subj_Level.id(1:13) = [];
for f = fieldnames(CM.wh_keep)'
    CM.wh_keep.(f{1})(1:13) = []; 
end


save(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')


%% Ns in each group
Ns = histc(get_var(CM, 'Group'), 1:3);
fprintf('Groups: %s\n', CM.Subj_Level.descrip{1})
fprintf('Group %i: N = %d\n', [[1:3]' Ns]')


%% robust regression WITHIN each group, using only Ss for that group. Basically a robust one-sample T-test
grp = get_var(CM, 'Group');
mygroup = 3; %edit to be 1, 2, or 3
within_group = get_wh_image(delta_relisten, grp==mygroup);
within_group.X = ones(Ns(mygroup),1); % estimate intercept (ttest)
out = regress(within_group, .001, 'unc', 'robust'); %  will do a robust one sample ttest within group, effectively
save(fullfile(outdir, ['within_group' num2str(mygroup) '_robust.mat']), 'out')


%% robust regression: CM vs. Faml and CM vs. OxyPla comparisons
dat = delta_relisten;

% build design matrix. Make CM the reference group (intercept), and add dummy regressors for Faml and OxyPla
dm = zeros(size(delta_relisten.dat, 2), 2);
dm( grp==2, 1) = 1;
dm( grp==3, 2) = 1;
figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm; % regs are, in order: Oxy vs CM, Faml vs CM, abs change in CM

%% estimate model
out2 = regress(dat, .001, 'unc', 'robust'); % automatically adds intercept as last column --> CM
save(fullfile(outdir, 'cm_vs_oxy_faml_robust.mat'), 'out2')


%% estimate placebo effects (OxyPla vs. Faml)
dat = delta_relisten;
grp = get_var(CM, 'Group');

% build design matrix. Make Faml the reference group (intercept), and add dummy regressors for CM and OxyPla
dm = zeros(size(delta_relisten.dat, 2), 2);
dm( grp==1, 1) = 1;
dm( grp==2, 2) = 1;
figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm; % regs are, in order: Cm vs Faml, Oxy vs Faml, abs change in Faml

%% estimate model
out2 = regress(dat, .001, 'unc', 'robust'); % automatically adds intercept as last column --> CM
save(fullfile(outdir, 'cm_oxy_vs_faml_robust.mat'), 'out2')









%% -- BELOW IS ALL OLD, TO BE DELETED



%% just to be double sure, estimate w/ contrast codes
dat = delta_relisten;

% contrast codes for CM vs. Faml, Oxy vs. other two
dm = zeros(size(delta_relisten.dat, 2), 2);
dm( grp==1, 1) = 1;
dm( grp==3, 1) = -1;
dm( grp==2, 2) = 1;
dm( grp~=2, 2) = -1;

figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm; % regs are, in order: Oxy vs CM, Faml vs CM, abs change in CM

% estimate model
out3 = regress(dat, .01, 'unc')%, 'robust'); % automatically adds intercept as last column --> CM
%save(fullfile(outdir, 'cm_vs_oxy_faml.mat'), 'out2')

%% check that CC and dummy give same results -- yes they do
a = get_wh_image(out3.t, 1); % CM vs. faml
b = get_wh_image(out2.t, 2); % Faml vs. CM
corr(a.dat, b.dat)
figure; scatter(a.dat, b.dat * -1);
diffs = abs(a.dat - b.dat*-1); max(diffs)

%% confirm that the within group estimate of CM effects, and the intercept from the 2nd
% model, give very similar results. will be slightly different b/c different number of Ss in both
% models

b = get_wh_image(out2.t,3); corr(out.t.dat, b.dat)
figure; scatter(out.t.dat, b.dat)


%% build model: group by time interaction, 
dat = delta_relisten;
grp = get_var(CM, 'Group');


dm = ones( size(delta_relisten.dat, 2), 1) * -1;
dm( grp==1, 1) = 1;
dm( grp==2, 1) = 0;
%dm = scale(dm, 1);
figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm;

%% estimate model
out = regress(dat, .001, 'unc');%, 'robust'); % automatically adds intercept as last column --> Faml
save(fullfile(outdir, 'CM_vs_Faml_norobust.mat'), 'out')


%% build model: group by time interaction, CM vs. combined controls

% build design matrix of 1 = CM, -1 = Faml and OxyPla
% intercept is avg change over time
dm = ones( size(delta_relisten.dat, 2), 1) * -1;
dm( grp==1, 1) = 1;
dm = scale(dm,1);
figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm;

%% estimate model
out = regress(dat, .001, 'unc', 'robust'); % automatically adds intercept as last column --> Faml
save(fullfile(outdir, 'CM_vs_combined_controls_robust.mat'), 'out')



%% build model: group by time interaction, CM vs. PLA

dat = delta_relisten;
grp = get_var(CM, 'Group');

% build design matrix of 1 = CM, -1 = OxyPla.
% intercept is avg change over time across all Ps, b/c I scale it. That's
% right -- the only impact of scaling will be on the interpretation of the
% intercept, whether it is avg effect over time OR effect in oxypla group.
% but scaling should have no effect on the comparison of interest. so can
% just just leave it in for consistency w/ the above
dm = zeros( size(delta_relisten.dat, 2), 1);
dm( grp==1, 1) = 1;
dm( grp==2, 1) = -1;
dm = scale(dm, 1);
figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm;

%% estimate model
out = regress(dat, .001, 'unc', 'robust'); % automatically adds intercept as last column --> Faml
save(fullfile(outdir, 'CM_vs_OxyPla_robust.mat'), 'out')




%% build model: group by time interaction, PLA vs. FAML

dat = delta_relisten;
grp = get_var(CM, 'Group');

% build design matrix of 1 = CM, -1 = OxyPla.
% intercept is avg change over time across all Ps, b/c I scale it. That's
% right -- the only impact of scaling will be on the interpretation of the
% intercept, whether it is avg effect over time OR effect in oxypla group.
% but scaling should have no effect on the comparison of interest. so can
% just just leave it in for consistency w/ the above
dm = zeros( size(delta_relisten.dat, 2), 1);
dm( grp==2, 1) = 1;
dm( grp==3, 1) = -1;
dm = scale(dm, 1);
figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm;

%% estimate model
out = regress(dat, .001, 'unc', 'robust'); % automatically adds intercept as last column --> Faml
save(fullfile(outdir, 'CM_vs_OxyPla_robust.mat'), 'out')










%% save the sgACC and parahip seeds for submission to connectivity analyses
load(fullfile(outdir, 'resultsCM_v_Faml_robust.mat'));
regions = region(out.t(1));
orthviews(regions(4)) % 4th region happens to be sgACC

sgACC = region2fmri_data(regions(4), out.t(1));
write(sgACC, 'fname', fullfile(outdir, 'sgACCmask_CMvsFamlrobust_001unc.nii'));

sgACC = region2fmri_data(regions(2), out.t(1));
write(sgACC, 'fname', fullfile(outdir, 'para-hip_mask_CMvsFamlrobust_001unc.nii'));



%% OLD: drop Ps from behavioral data who are excluded.  run this cell once.

% NO -- don't drop. No empirical reason to think they are outliers (mahal and image mean all ?normal?), so include them.
% drop P who did not understand donation instructions (29) and P with a
% technical error with biography timings (17)
%to_drop = [17 29];
%CM.Subj_Level.data(to_drop, :) = [];
%CM.Subj_Level.id(to_drop) = [];

% must drop them from brain data too. ALREADY DONE in the data folders
%wh = true(size(time1.dat, 2), 1);
%wh(to_drop) = 0;
%time1 = get_wh_image(time1, wh);
%time2 = get_wh_image(time2, wh);
