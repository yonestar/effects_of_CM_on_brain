clear all

%% load data and other set up
 
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
outdir = fullfile(repodir, 'results', 'group_level');
cd(repodir);

doplot=1;

time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

gm_mask = fmri_data(which('gray_matter_mask.img'));
%gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));
gm_mask.dat = gm_mask.dat > .25; % tighter gm mask. increase val to make tighter

time1 = apply_mask(time1, gm_mask);
time2 = apply_mask(time2, gm_mask);

load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')


% drop P who misunderstood donation instructions
wh = get_var(CM, 'deltadon') > -45 | get_var(CM, 'deltadon') < -46 ;

time1 = get_wh_image(time1, wh);
time2 = get_wh_image(time2, wh);
delta_relisten = image_math(time2, time1, 'minus');

grp = get_var(CM, 'Group',wh);
grp_cc(:,1) = (grp == 1) + -1*(grp ~= 1);
grp_cc(:,2) = (grp == 2) + -1*(grp == 3);

% get deltadon/deltaFAS and mean center within group
deltadon_mcgrp = get_var(CM,'deltadon', wh);
don1_mcgrp = get_var(CM,'don1', wh);
don2_mcgrp = get_var(CM,'don2', wh);
deltaFAS_mcgrp = get_var(CM,'deltaFAS', wh);

for i=1:3
    deltadon_mcgrp(grp==i) = scale(deltadon_mcgrp(grp==i), 1);
    don1_mcgrp(grp==i) = scale(don1_mcgrp(grp==i), 1);
    don2_mcgrp(grp==i) = scale(don2_mcgrp(grp==i), 1);
    deltaFAS_mcgrp(grp==i) = scale(deltaFAS_mcgrp(grp==i), 1);
end

% create a don2 orthogonalized wrt to don1
mdl = fitglm(don1_mcgrp, don2_mcgrp);
don2_mcgrp_orth = mdl.Residuals.Raw;

if doplot
    figure;
    obj = canlab_results_fmridisplay('compact2', 'overlay', 'icbm152_2009_symmetric_for_underlay.img' );
end

%% robust regression: CM vs. Faml and CM vs. OxyPla comparisons, with DeltaDon and DeltaFAS mean centered within group
% deltadon and deltaFAS only share 3.5% of variance after mean-centering
% within group, so OK include in same model. Also with robust corr, not a
% strong relationship between the two vars (

dat = delta_relisten;

% build design matrix. 
% Make CM the reference group (intercept), and add dummy regressors for Faml and OxyPla
dm = zeros(size(delta_relisten.dat, 2), 2);
dm( grp==2, 1) = -1;
dm( grp==3, 2) = -1;

% add deltadon/fas mean centered within grp, or don1/2
dm(:,3) = deltadon_mcgrp;
%dm(:,3) = don1_mcgrp;
%dm(:,4) = don2_mcgrp_orth;

%figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm; % regs are, in order: Oxy vs CM, Faml vs CM, abs change in CM

% estimate model
out2 = regress(dat, .001, 'unc', 'k', 10); % automatically adds intercept as last column --> CM
%save(fullfile(outdir, 'cm_vs_oxy_faml_deltadonfas_robust.mat'), 'out2')



%% which region survive for deltadon?  L NAc, pOFC, others

if doplot
    print_header('Delta donation predicting delta brain');
    toshow = get_wh_image(out2.t,3);
    obj = removeblobs(obj);
    obj = addblobs(obj, region(toshow));
    snapnow
end

table(region(get_wh_image(out2.t, 3)));


%% So, no regions positively relate to delta don cntrolling for group, while a number of mOFC/NAc regions *negatively* relate to delta don. Let us explore this further by decomposing "delta brain" and deltadon into Time1 and Time2 variables.


%% Predict Time1 brain from Don1 + Don2 orthogonalized wrt to Don1 + group

dat = time1;

% build design matrix. 
% Make CM the reference group (intercept), and add dummy regressors for Faml and OxyPla
dm = zeros(size(delta_relisten.dat, 2), 2);
dm( grp==2, 1) = -1;
dm( grp==3, 2) = -1;

% add deltadon/fas mean centered within grp, or don1/2
dm(:,3) = don1_mcgrp;
dm(:,4) = don2_mcgrp_orth;

%figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm; % regs are, in order: Oxy vs CM, Faml vs CM, abs change in CM

% estimate model
out2 = regress(dat, .001, 'unc', 'k', 10); % automatically adds intercept as last column --> CM
%save(fullfile(outdir, 'cm_vs_oxy_faml_deltadonfas_robust.mat'), 'out2')


%% Strongest results are for Don2 (orthogonalized wrt to Don1) predicting brain response at Time1
if doplot
    print_header('Don2 (orthogonalized wrt to Don1) predicting Time 1 brain');
    toshow = get_wh_image(out2.t,4);
    obj = removeblobs(obj);
    obj = addblobs(obj, region(toshow));
    snapnow
end

table(region(get_wh_image(out2.t, 4)));

%% Don1 is not strongly related to any brain responses at Time1, though
if doplot
    print_header('Don1 predicting Time 1 brain');
    toshow = get_wh_image(out2.t,3);
    obj = removeblobs(obj);
    obj = addblobs(obj, region(toshow));
    snapnow
end

table(region(get_wh_image(out2.t, 3)));

%% Brain activity at Time1 in mOFC/NAc regions prospectively predicts who will increase or decrease in donation, across groups. Like a prognostic marker of natural history of donation change

%% Is this prognostic finding driven by a particular group? Test this w/ a group * Don2 interaction term in the model
% Interactions terms are between the two group contrast codes (CM vs. cntrls, OxyPla vs. Faml) and Don2 orthogonalized wrt to Don1. 

% build design matrix. 
% Make CM the reference group (intercept), and add dummy regressors for Faml and OxyPla
dm = grp_cc;

% add deltadon/fas mean centered within grp, or don1/2
dm(:,3) = don1_mcgrp;
dm(:,4) = don2_mcgrp_orth;

dm(:,5) = don2_mcgrp_orth .* dm(:,1);
dm(:,6) = don2_mcgrp_orth .* dm(:,2);

%figure; imagesc([dm ones(length(dm), 1)]); colorbar
dat.X = dm; % regs are, in order: Oxy vs CM, Faml vs CM, abs change in CM

% estimate model
out2 = regress(dat, .001, 'unc', 'k', 10); % automatically adds intercept as last column --> CM
%save(fullfile(outdir, 'cm_vs_oxy_faml_deltadonfas_robust.mat'), 'out2')

%% Results: When these interaction terms are included in the model, much less survives threshold. This suggests that mOFC/NAC predicts Don2 across groups