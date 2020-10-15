clear all

%% load data and other set up
 
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
outdir = fullfile(repodir, 'results', 'group_level');
rdir = fullfile(repodir, 'data', 'for_randomise');

time1 = fmri_data(filenames(fullfile(repodir, 'data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile(repodir, 'data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

% time2 - time1
delta_relisten = image_math(time2, time1, 'minus'); 

%% load group assignments and behavioral data from trial
% this data is sorted by subject num in alphanumeric order
% so it maps on to the imaging data one-to-one, in order

load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')
grp = get_var(CM, 'Group');

%% make masks for randomise
atl = load_atlas('canlab2018_2mm');
subcortical = load_atlas('cit168');

% BA 10 = mPFC
mkdir(fullfile(rdir, 'masks', 'mPFC'));
mPFC = select_atlas_subset(atl, {'Ctx_10'}, 'flatten');

mkdir(fullfile(rdir, 'masks', 'ofc'));
ofc = select_atlas_subset(atl, {'OFC'}, 'flatten');
write(ofc, 'fname', fullfile(rdir, 'masks', 'ofc', 'ofc.nii'), 'mni', 'overwrite')

%%
%combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(subcortical, {'VTA'}));
mkdir(fullfile(rdir, 'masks', 'nac'));
nac =  select_atlas_subset(subcortical, {'NAC'});

% resample nac

% resample NAC from CIT168 to 2x2x2, without left-shifting
% flirt -in nac.nii -ref nac.nii -out nac_RS -applyisoxfm 2
nac_rs = fmri_data('nac_RS_MNI_flirt.nii.gz')
% orthviews_multiple_objs({nac, nac_rs})

% its still a bit shifted, but whatever -- i tried!!
write(nac_rs, 'fname', fullfile(rdir, 'masks', 'nac', 'nac_rs_MNI.nii'), 'overwrite','mni');

%% combined NAC / OFC

nac = fmri_data(fullfile(rdir, 'masks', 'nac', 'nac_rs_MNI.nii'));
ofc = fmri_data(fullfile(rdir, 'masks', 'ofc', 'ofc.nii'));
nac_ofc = image_math(nac, ofc, 'plus');
nac_ofc.dat = nac_ofc.dat>0;
write(nac_ofc, 'fname', fullfile(rdir, 'masks', 'nac_ofc', 'nac_ofc.nii'), 'overwrite','mni');


%% write out dat_delta, design, and contrasts -- will be shared for all masks!

% write(delta_relisten, 'fname', fullfile(rdir, 'delta_relisten.nii'), 'mni', 'overwrite');

delta_relisten_CMvPLA = get_wh_image(delta_relisten, grp ~= 3);
write(delta_relisten_CMvPLA, 'fname', fullfile(rdir, 'delta_relisten_CMvPLA.nii'), 'mni', 'overwrite');

delta_relisten_CMvFAML = get_wh_image(delta_relisten, grp ~= 2);
write(delta_relisten_CMvFAML, 'fname', fullfile(rdir, 'delta_relisten_CMvFAML.nii'), 'mni', 'overwrite');

%% Design: CM vs PLA
grp_CMvPLA = grp(grp~=3);
X = double(grp_CMvPLA==1);
X(X==0) = -1;

writetable(array2table(X), fullfile(rdir, 'design_CMvPLA.txt'), 'WriteVariableNames', false, 'delimiter', ' ');
cd(rdir)
!/usr/local/fsl/bin/Text2Vest design_CMvPLA.txt design_CMvPLA.mat

%% Design: CM vs FAML
grp_CMvFAML = grp(grp~=2);
X = double(grp_CMvFAML==1);
X(X==0) = -1;

writetable(array2table(X), fullfile(rdir, 'design_CMvFAML.txt'), 'WriteVariableNames', false, 'delimiter', ' ');
cd(rdir)
!/usr/local/fsl/bin/Text2Vest design_CMvFAML.txt design_CMvFAML.mat


%% contrast: will always be the same (CM > control). One-sided test -- OK.
f = fopen(fullfile(rdir, 'contrast.txt'), 'w');
fprintf(f, '1\n');
fclose(f);
!/usr/local/fsl/bin/Text2Vest contrast.txt contrast.con
