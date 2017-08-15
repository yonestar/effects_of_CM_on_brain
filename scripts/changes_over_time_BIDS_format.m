%% load data
% data is in BIDS format (so can later share on openfmri)
% right now, only have derivatives (the pre and post contrast images) in this format
% later, put the full data in BIDS format for sharing on openfMRI
 
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
dm(wh_oxy, 2) = 1;

dat.X = dm;

dat = apply_mask(dat, gm_mask);

out = regress(dat, .001, 'unc')%, 'robust'); % automatically adds intercept as last column --> Faml

%% save output files
outdir = fullfile(basedir, 'derivatives', 'group level');

save(fullfile(outdir, 'resultsEachGroup_robust.mat'), 'out')

write(get_wh_image(out.t,1), 'fname', 'deltarelistenCM_robust.nii', 'mni', 'thresh')

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
dat = dats.listen_delta;

dm = get_var(CM, 'deltaFAS', CM.wh_keep.t2);

dm(4) = 0; % this P is to be excluded.

% mean center within group
dm(wh_med == 1) = detrend(dm(wh_med == 1), 'const');
dm(wh_oxy == 1) = detrend(dm(wh_oxy == 1), 'const');
dm(wh_faml == 1) = detrend(dm(wh_faml == 1), 'const');

dat.X = dm;
out=regress(dat, .005, 'unc');

%multi_threshold(get_wh_image(out.t,1), 'thresh', [.005 .01 .05]) 



%% differences at time2, controlling for time1
datT2 = dats.relistenT2;
datT2 = apply_mask(datT2, which('gray_matter_mask.img'));

datT1 = apply_mask(dats.relistenT1, which('gray_matter_mask.img'));

% intercept for each group
dm = zeros(57,3);
dm(wh_med, 1) = 1;
dm(wh_oxy, 2) = 1;
dm(:, 3) = 1; %faml

% CM vs. Oxy
dm = zeros(57,2);
dm(wh_med, 1) = 1;
dm(wh_oxy, 1) = -1;
dm(:, 2) = 1; %faml

% Oxy vs. Faml
dm = zeros(57,2);
dm(wh_oxy, 1) = 1;
dm(wh_faml, 1) = -1;
dm(:, 2) = 1; %CM


n = 57;

%%
for i = 1:size(datT2.dat,1)
        
    % Create X    
    X = [dm datT1.dat(i, :)'];
    Y = datT2.dat(i,:)';
    
  %  [n, k] = size(X);
    [bb,stats] = robustfit(X, Y, 'bisquare', [], 'off');

    b(:,i)=bb; %Betas
    t(:,i)=stats.t; %t-values
    p(:,i)=stats.p; %p-values
    dfe(:,i)=stats.dfe; %degrees of freedom
    stderr(:,i)=stats.se; %Standard Error
    sigma(:,i)=stats.robust_s; %robust estimate of sigma. LC not sure this is the best choice can switch to 'OLS_s','MAD_s', or 's'
    %r(:,i) = dat.Y - X * b(:,i); %residual
end

% save results
out = struct;
inputargs = {.01, 'unc'};
dat = datT1;

% Betas
out.b = statistic_image;
out.b.type = 'Beta';
out.b.p = p';
out.b.ste = stderr';
out.b.N = n;
out.b.dat = b';
out.b.dat_descrip = sprintf('Beta Values from regression, intercept is last');
out.b.volInfo = dat.volInfo;
out.b.removed_voxels = dat.removed_voxels;
out.b.removed_images = false;  % this image does not have the same dims as the original dataset
out.b = threshold(out.b, inputargs{:}); % Threshold image

% T stats
out.t = statistic_image;
out.t.type = 'T';
out.t.p = p';
out.t.ste = stderr';
out.t.N = n;
out.t.dat = t';
out.t.dat_descrip = sprintf('t-values from regression, intercept is last');
out.t.volInfo = dat.volInfo;
out.t.removed_voxels = dat.removed_voxels;
out.t.removed_images = false;  % this image does not have the same dims as the original dataset
out.t = threshold(out.t, inputargs{:}, 'noverbose'); %Threshold image

% DF as fmri_data
out.df = dat;
out.df.dat = dfe';
out.df.dat_descrip = sprintf('Degrees of Freedom');

% Sigma as fmri_data
out.sigma = dat;
out.sigma.dat = sigma';
out.sigma.dat_descrip = sprintf('Sigma from Regression');


save(fullfile(basedir, 'results.mat'), 'out')

%% view montages with multi threshold

multi_threshold(get_wh_image(out.t,1), 'thresh', [.005 .01 .05]) 


%% view on full montage

ofull = canlab_results_fmridisplay('montagetype', 'full');
%%

dat = get_wh_image(out.t, 2);

dat = threshold(dat, .005, 'unc');
red = [237,95,52]/256; 
ofull = addblobs(ofull, region(dat) , 'splitcolor', {[0 0 1] [.3 0 .8] red-.2 red+.2}, 'wh_montages', 1);

dat = threshold(dat, .001, 'unc');
ofull = addblobs(ofull, region(dat) , 'splitcolor', {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]}, 'wh_montages', 1);

%%
ofull = removeblobs(ofull)


%% RM-ANOVA.  this is like the pdf phil sent, but don't see other refs for it.

% concatenate all the data
% add a regressor for group
% a regressor for time (1 or 2)
% an intercept per subject
% a group by time interaction
% 

%% non-parametric (in vmPFC)
%R = robust_reg_nonparam(dm, 5, 'mask', which('VMPFC_right.img'), 'names', {'cm_vs_faml', 'intcp'}, 'data', dat)


%% try my own SVC version.  FDR in vmPFC
dat = apply_mask(dat, vmpfc);
%orthviews(get_wh_image(dat,1))
dat.X = dm;
out=regress(dat, .05, 'fdr');
