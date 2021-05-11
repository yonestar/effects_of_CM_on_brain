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


%% make mesolimbic masks
% cons: NAC is shifted. does not include pgACC areas
% pros: VTA is fine. Includes the regions in my original submission

% VTA: when I load it from 'subcortical' atlas, it gets pushed around badly in the resampling. But when I load it from 'canlab2018_2mm', its much larger but centered
% however, this doesn't help for NAC -- its still shifted. Don't know how
% to get a not-shifted NAC in 2mm space.

atl = load_atlas('canlab2018_2mm');
combined_atlas = select_atlas_subset(atl, {'Ctx_10', 'OFC', 'NAC', 'VTA'}, 'flatten');
write(combined_atlas, 'fname', fullfile(rdir, 'masks', 'mesolimbic', 'mesolimbic.nii'), 'overwrite','mni'); % NAC is shifted but this may be the best I can do
mesolimbic = fmri_data(fullfile(rdir, 'masks', 'mesolimbic', 'mesolimbic.nii'));

combined_atlas = select_atlas_subset(atl, { 'OFC', 'NAC', 'VTA'}, 'flatten');
write(combined_atlas, 'fname', fullfile(rdir, 'masks', 'ofc_nac_vta', 'ofc_nac_vta.nii'), 'overwrite','mni'); % NAC is shifted but this may be the best I can do

%% make dmPFC / TPJ masks

% Mascaro ROIs: dmPFC [-9 50 37] mentalizing; IFG: [-45 29 4] mirror or control
% Weng ROIs:   IPC/TPJ: [46, -62, 36] mentalizing or mirror;   dlPFC: [38, 22, 46] mirror or control

mm = [-9 50 37; 9 50 37];
indx = iimg_xyz2spheres(mm2voxel(mm,time1.volInfo),time1.volInfo.xyzlist,5); % make 10mm radius sphere. Our voxels here at 2x2x2, so 3 vox
dmPFC = get_wh_image(time1, 1); % arbitrary
dmPFC.dat = indx;

mm = [46, -62, 36; -46, -62, 36];
indx = iimg_xyz2spheres(mm2voxel(mm,time1.volInfo),time1.volInfo.xyzlist,5); % make 10mm radius sphere. Our voxels here at 2x2x2, so 3 vox
TPJ = get_wh_image(time1, 1); % arbitrary
TPJ.dat = indx;

dmpfc_tpj = image_math(dmPFC, TPJ, 'plus');
orthviews(dmpfc_tpj)

write(dmpfc_tpj, 'fname', fullfile(rdir, 'masks', 'dmpfc_tpj', 'dmpfc_tpj.nii'), 'overwrite','mni');


% Could take take it from an atlas, but its 10x the size of ofc_nac_vta

 % confirmed that DefaultA includes Weng coords
 % confirmed that DefaultB includes Mascaro coords



%% write out two dat_deltas

delta_relisten_CMvPLA = get_wh_image(delta_relisten, grp ~= 3);
write(delta_relisten_CMvPLA, 'fname', fullfile(rdir, 'delta_relisten_CMvPLA.nii'), 'mni', 'overwrite');

delta_relisten_CMvFAML = get_wh_image(delta_relisten, grp ~= 2);
write(delta_relisten_CMvFAML, 'fname', fullfile(rdir, 'delta_relisten_CMvFAML.nii'), 'mni', 'overwrite');

%% Write Design: CM vs PLA
grp_CMvPLA = grp(grp~=3);
X = double(grp_CMvPLA==1);
X(X==0) = -1;

writetable(array2table(X), fullfile(rdir, 'design_CMvPLA.txt'), 'WriteVariableNames', false, 'delimiter', ' ');
cd(rdir)
!/usr/local/fsl/bin/Text2Vest design_CMvPLA.txt design_CMvPLA.mat

%% Write Design: CM vs FAML
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



%% Randomize command:
% randomise -i ../../../delta_relisten_CMvFAML.nii -o CMvFAML -d ../../../design_CMvFAML.mat -t ../../../contrast.con -m ../ofc_nac_vta.nii -n 5000 -D -T -v 6
% randomise -i ../../../delta_relisten_CMvPLA.nii -o CMvPLA -d ../../../design_CMvPLA.mat -t ../../../contrast.con -m ../ofc_nac_vta.nii -n 5000 -D -T -v 6

% to view Range of (1-p) values:
% fslstats <corr p img> -R

% to create clusters image:
% cluster -i CMvFAML_tfce_corrp_tstat1.nii.gz -t .95 -c CMvFAML_tstat1.nii.gz --scalarname="1-p" --mm -o clusters_img > clusters_corrp.txt


%% view output at .05
% p=statistic_image('grpbytime_tfce_corrp_tstat1.nii', 'type', 'p');
% cl = fmri_data('clusters_img.nii.gz')
% orthviews_multiple_objs({threshold(p, [.95 Inf], 'raw-between'), cl})
