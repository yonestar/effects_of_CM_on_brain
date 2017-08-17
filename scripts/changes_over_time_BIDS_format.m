%% load data
% data is in BIDS format (so can later share on openfmri)
% right now, only have derivatives (the pre and post contrast images) in this format
% later, put the full data in BIDS format for sharing on openfMRI
 
repo_dir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
outdir = fullfile(basedir, 'derivatives', 'group level');
cd(repo_dir);

% load group assignments and behavioral data from trial
load CM_data_2015-07-29_16_16.mat  % canlab dataset object.  loads 'CM' var

% load brain data.  point to contrast images on local computer (for now)
basedir = '/Users/yoni/Dropbox/Research/LKM/projects/LKM_fmri_longitudinal/Imaging/Analyses/BIDS';
cd(fullfile(basedir, 'derivatives', 'subject level'))

% con_0011 is the change in [relisten > arrows] contrast, pre-to-post
delta_relisten_fnames = filenames(fullfile('*', 'con_0011.img'));
delta_relisten = fmri_data(delta_relisten_fnames);

% gray matter mask
gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));


%% make wh_keeps as they should be for these analyses

wh_med = CM.wh_keep.med;
wh_oxy = CM.wh_keep.oxy;
wh_faml = CM.wh_keep.faml;

% drop the P with bad bios timing
wh_faml( find(~CM.wh_keep.bios_timing_correct) ) = 0;

% add back in the P I dropped b/c of wierd donation behavior
%ind = find( xor(CM.wh_keep.t2, CM.wh_keep.t2dropone) );
%wh_med(ind) = 1;

% now drop them down to a space of 57 Ps (the t2 space)
wh_med = wh_med(CM.wh_keep.t2);
wh_oxy = wh_oxy(CM.wh_keep.t2);
wh_faml = wh_faml(CM.wh_keep.t2);

% check for correctness
sum(wh_med), sum(wh_oxy), sum(wh_faml)
groups = double([wh_med wh_oxy wh_faml]);


%% within group change over time
dat = delta_relisten;

dm = zeros(57,2);
dm(wh_med, 1) = 1;
dm(wh_faml, 2) = -1;
dm(wh_oxy, 2) = 0;

dat.X = dm;

dat = apply_mask(dat, gm_mask);

out = regress(dat, .005, 'unc')%, 'robust'); % automatically adds intercept as last column --> Faml

%% save output files

save(fullfile(outdir, 'resultsEachGroup_robust.mat'), 'out')

write(get_wh_image(out.t,1), 'fname', 'deltarelistenCM_robust.nii', 'mni', 'thresh')

% save output as masks for use as seed in connectivity analyses

%% robust regression for the 3 group x time interactions of interest

dat = apply_mask(dats.relisten_delta, gm_mask);

% CM v Oxy
dm = zeros(57,1);
dm(wh_med) = 1 / sum(wh_med);
dm(wh_oxy) = -1 / sum(wh_oxy);
dat.X = dm;

out = regress(dat, 'robust'); % automatically adds intercept as last column
save(fullfile(outdir, 'resultsCM_v_Oxy_robust.mat'), 'out', 'dat') % save output file

% CM v Faml
dm = zeros(57,1);
dm(wh_med) = 1 / sum(wh_med);
dm(wh_faml) = -1 / sum(wh_faml);
dat.X = dm;

out = regress(dat, 'robust'); % automatically adds intercept as last column
save(fullfile(outdir, 'resultsCM_v_Faml_robust.mat'), 'out', 'dat') % save output file

% Oxy v Faml
dm = zeros(57,1);
dm(wh_oxy) = 1 / sum(wh_oxy);
dm(wh_faml) = -1 / sum(wh_faml);
dat.X = dm;

out = regress(dat, 'robust'); % automatically adds intercept as last column
save(fullfile(outdir, 'resultsOxy_v_Faml_robust.mat'), 'out', 'dat') % save output file

%% view CM vs. Faml effect at multiple thresholds

load(fullfile(outdir, 'resultsCM_v_Faml_robust.mat'))

[~, sig] = multi_threshold(get_wh_image(out.t,1), 'thresh', [.001 .005 .01], 'sizethresh', [1 1 1], 'nodisplay');
% returns *sig:**
%        vector of significant voxels at each thresh, for each region
%               - cell array of images in object with matrix of values
%               for each threshold

% save three different output files (at each thresh) for mricrogl
a = replace_empty(get_wh_image(out.t,1));

a.sig = sig{1}(:,1);
write(a, 'mni', 'thresh', 'fname', fullfile(outdir, 'CMvFaml_001unc.nii'));
a.sig = sig{1}(:,2);
write(a, 'mni', 'thresh', 'fname', fullfile(outdir, 'CMvFaml_005unc.nii'));
a.sig = sig{1}(:,3);
write(a, 'mni', 'thresh', 'fname', fullfile(outdir, 'CMvFaml_01unc.nii'));


%% save the sgACC seed for submission to connectivity analyses
load(fullfile(outdir, 'resultsCM_v_Faml_robust.mat'));
regions = region(out.t(1));
orthviews(regions(4)) % 4th region happens to be sgACC

sgACC = region2fmri_data(regions(4), out.t(1));
write(sgACC, 'fname', fullfile(outdir, 'sgACCmask_CMvsFamlrobust_001unc.nii'));

%% print cluster table of CM vs. Faml results
cluster_table(out.t(1)) % fails w/ error









%% SEMI-OLD   display results of main effects and interactions

%o2 = canlab_results_fmridisplay('montagetype', 'compact2');

o2 = montage(o2, 'axial', 'wh_slice', [-15 -15 -15], 'noverbose');

% shift all axes down and right
allaxh = o2.montage{1}.axis_handles;
for i = 1:length(allaxh)
    pos1 = get(allaxh(i), 'Position');
    pos1(2) = pos1(2) - 0.08;
    pos1(1) = pos1(1) + 0.03;
    set(allaxh(i), 'Position', pos1);
end

axh = axes('Position', [0.0 0.08 .15 1]);
o2 = montage(o2, 'saggital', 'wh_slice', [8 8 8], 'spacing', 2, 'existing_axes', axh, 'noverbose');
enlarge_axes(gcf, 1.3);

o2 = multi_threshold(mydat, 'thresh', [.001 .005 .01], 'o2', o2, 'sizethresh', [1 1 1]);


% show the conjunction of the interaction and the absolute increase in CM group.  p
% values shown in figure should be for interaction.  so the "mask" is for
% the absolute increase, and is also pruned at [.001 .005 .01]. 

%% show effect of CM, and CM vs. faml interaction
o2 = removeblobs(o2);

% make mask
load resultsEachGroup_robust

[o2, sig] = multi_threshold(get_wh_image(out.t,1), 'thresh', [.001 .005 .01], 'o2', o2, 'sizethresh', [1 1 1], 'nodisplay');
mask = replace_empty(get_wh_image(out.t,1));
mask.sig = sig{1}(:,3);
mask.dat(mask.sig==0) = 0;
mask = remove_empty(mask);
orthviews(mask)

% mask interaction by within-group effect, and then display interaction
load resultsCM_v_Faml_robust
mydat = get_wh_image(out.t,1);
mydat = apply_mask(mydat, mask);

[o2, sig] = multi_threshold(mydat, 'thresh', [.001 .005 .01], 'o2', o2, 'sizethresh', [1 1 1], 'nodisplay');

% save three different output files (at each thresh) for mricrogl
a = replace_empty(mydat);

a.sig = sig{1}(:,1);
write(a, 'mni', 'thresh', 'fname', 'CMvFaml_001unc_masked_CM_abs_increase.nii');
a.sig = sig{1}(:,2);
write(a, 'mni', 'thresh', 'fname', 'CMvFaml_005unc_masked_CM_abs_increase.nii');
a.sig = sig{1}(:,3);
write(a, 'mni', 'thresh', 'fname', 'CMvFaml_01unc_masked_CM_abs_increase.nii');







%% OLD show an outline
outlinevars = {'outline', 'linewidth', 3, 'trans', 'transvalue', 0, 'color', [1 1 1]};
o2 = addblobs(o2, region(threshold(get_wh_image(out.t, 1), .005, 'unc')), outlinevars{:});


%% within-group changes in FAS/deltadon

dm = get_var(CM, 'deltaFAS', CM.wh_keep.t2);

dm(4) = 0; % this P is to be excluded.

% mean center within group
dm(wh_med == 1) = detrend(dm(wh_med == 1), 'const');
dm(wh_oxy == 1) = detrend(dm(wh_oxy == 1), 'const');
dm(wh_faml == 1) = detrend(dm(wh_faml == 1), 'const');

dat.X = dm;
out=regress(dat, .005, 'unc');

%multi_threshold(get_wh_image(out.t,1), 'thresh', [.005 .01 .05]) 

