%{
my summary to self: 

Final model I decided on was:
 T2 ROI value = T1 ROI value + deltadon mean-centered within group (i.e., controlling for group) + intcpt
 
I ran this for all 6 ROIs. I also ran 6 models with deltaFAS (instead of
deltadon).

For all ROIs, there was no effect of deltadon or deltaFAS on the T2 ROI value when
controlling for group, with one exception across those 12 tests.

The exception: T2 NAc was *negatively* predicted by deltadon. 
To explore this further, I ran this model:
 
T2 NAC = T1 NAC + T1 donation (mean-centered within grp) + T2 donation (mean-centered within grp)
 
Found that T1 donation positively predicted T2 NAc, while T2 donation negatively
predicted T2 NAc. This is odd, and, explains why deltadon negatively
predicts T2 NAc.
%}

clear all

%%
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
figdir = fullfile(repodir, 'figures');
outdir = fullfile(repodir, 'results', 'group_level');
cd(repodir);


load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')

wh = get_var(CM, 'deltadon') > -45 | get_var(CM, 'deltadon') < -46 ; % drop the P who misunderstood donation instructions

grp = get_var(CM, 'Group', wh);
don1 = get_var(CM, 'don1', wh);
don2 = get_var(CM, 'don2', wh);
fas1 = get_var(CM, 'FAS1', wh);
fas2 = get_var(CM, 'FAS2', wh);

deltadon = don2 - don1;
deltaFAS = fas2 - fas1;

time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

% drop P who misunderstood donation instructions
time1 = get_wh_image(time1, wh);
time2 = get_wh_image(time2, wh);
delta_relisten = image_math(time2, time1, 'minus');

colors = {[1 .2 .2] [.2 .2 1] [.4 .4 .4] [.7 .7 .7]};

grp_cc(:,1) = (grp == 1) + -1*(grp ~= 1);
grp_cc(:,2) = (grp == 2) + -1*(grp == 3);


%% make my atlas with regions I want
subcortical = load_atlas('subcortical_rl');
glasser = load_atlas('glasser');

% BA 10 = mPFC
combined_atlas = merge_atlases(select_atlas_subset(glasser, {'Ctx_10'}, 'flatten'), select_atlas_subset(glasser, {'OFC'}, 'flatten'));
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(subcortical, {'NAC'}));
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(subcortical, {'VTA'}));

% dmPFC (9m) and TPJ (PFm
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(glasser, {'9m'}, 'flatten'));
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(glasser, {'PFm'}, 'flatten'));

orthviews(combined_atlas);
% mPFC{i} = vals(1).dat;
% OFC{i} = vals(2).dat;
% nac{i} = vals(3).dat;
% vta{i} = vals(4).dat;
% 
% dmPFC{i} = vals(5).dat;
% TPJ{i} = vals(6).dat;

%% extract from all ROIs for the 3 groups
roivals_T1 = extract_roi_averages(time1, combined_atlas);
roivals_T2 = extract_roi_averages(time2, combined_atlas);
roivals_delta = extract_roi_averages(delta_relisten, combined_atlas);

%% look at the data. 
create_figure('roi');
whroi = 3; % 3 = NAc; 2 = OFC -- the only regions showing a sig grp effect
scatterhist(roivals_T1(whroi).dat, roivals_T2(whroi).dat), lsline

create_figure('don');
scatterhist(don1, don2), lsline

create_figure('fas');
scatterhist(fas1, fas2), lsline

%% change in brain by change in behavior
clc
names = {'Intercept' 'deltaROI' 'CMvsCntrl' 'OxyVsFaml'};
whroi = 2; % 3 is NAc, 2 is mOFC

fitglm([roivals_delta(whroi).dat grp_cc], deltaFAS, 'VarNames', {names{2:end} 'deltaBehav'}) % non-robust -- see effects!

[~,stats]=robustfit([roivals_delta(whroi).dat grp_cc], deltaFAS); 

for i=1:length(names)
    fprintf('%s\tT = %3.2f\t p = %3.2f\n', names{i}, stats.t(i), stats.p(i));
end

%% full model at both timepoints
% predict time2 don/fas from baseline don/fas and pre and post ROI and grp

whroi = 2; % 3 is NAc, 2 is mOFC
[~,stats]=robustfit([don1 roivals_T1(whroi).dat roivals_T2(whroi).dat grp_cc], don2); 

clc
names = {'Intercept' 'T1 behav' 'T1 ROI  ' 'T2 ROI  ' 'CMvsCntrl' 'OxyVsFaml'};
for i=1:length(names)
    fprintf('%s\tT = %3.2f\t p = %3.2f\n', names{i}, stats.t(i), stats.p(i));
end

% time2 donation is predicted by higher NAC at time1, lower NAC at T2, more
% donation at T1. NO, its all outliers, all those effects disappear w/
% robust regression!!



%% WITHIN GROUP: T2 don pred by T1 don, change in NAC

i=1 % grp

% change in brain by change in behavior
clc
names = {'Intercept' 'deltaROI'};
whroi = 3; % 3 is NAc, 2 is mOFC

%fitglm([roivals_delta(whroi).dat grp_cc], deltaFAS, 'VarNames', {names{2:end} 'deltaBehav'}) % non-robust -- see effects!

[~,stats]=robustfit([roivals_delta(whroi).dat(grp==i)], deltadon(grp==i)); 

for i=1:length(names)
    fprintf('%s\tT = %3.2f\t p = %3.2f\n', names{i}, stats.t(i), stats.p(i));
end




%% WITHIN GROUP: plot ROI changes by change in don/FAS, within group

% choose var
whroi = 3; 

mdl = fitglm([roivals_delta(whroi).dat grp_cc], deltadon)

create_figure('scatter grp'); 
for i=1:3
    scatter(roivals_delta(whroi).dat(grp==i), deltadon(grp==i)), hold on
    h=lsline;
    set(h,'Color',colors{i})
    
end
xlabel('Change in ROI'), ylabel('Change in don')

