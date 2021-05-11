clear all

%% set up and load data
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
figdir = fullfile(repodir, 'figures');
gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));
underlay=which('keuken_2014_enhanced_for_underlay.img');
%overlay=which('SPM8_colin27T1_seg.img');
%overlay = which('icbm152_2009_symm_enhanced_for_underlay.img')

outlinecolor = [0 0 0];
splitcolor = {[0 0 1] [.3 0 .8] [.8 .3 0] [1 1 0]};

% load between group contrasts
load(fullfile(repodir,'results','group_level','cm_vs_oxy_faml_robust.mat'));
cm_vs_oxy = get_wh_image(out2.t, 1);
cm_vs_oxy.dat = cm_vs_oxy.dat * -1; % flip so it is cm vs oxy, instead of oxy vs cm
cm_vs_oxy = apply_mask(cm_vs_oxy, gm_mask); % apply GM mask

cm_vs_faml = get_wh_image(out2.t, 2);
cm_vs_faml.dat = cm_vs_faml.dat * -1; % flip so it is cm vs faml, instead of faml vs cm
cm_vs_faml = apply_mask(cm_vs_faml, gm_mask); % apply GM mask

%% for placebo effects
load(fullfile(repodir,'results','group_level','cm_oxy_vs_faml_robust.mat'));
oxy_vs_faml = get_wh_image(out2.t, 2);
oxy_vs_faml = apply_mask(oxy_vs_faml, gm_mask); % apply GM mask

%% NOTE: removeblobs() is not working correctly for surfaces, so need to recreate the montage in all the below
% another note: save surface brains as bitmaps; no advantage to svg (and,
% those white lines show up when saved as svg)

%% save out image files (e.g., for upload to NeuroVault)

write(cm_vs_faml, 'fname', fullfile(repodir,'results','group_level','cm_vs_faml.nii'))
write(cm_vs_oxy, 'fname', fullfile(repodir,'results','group_level','cm_vs_oxy.nii'))
write(oxy_vs_faml, 'fname', fullfile(repodir,'results','group_level','oxy_vs_faml.nii'))

i=1; load(fullfile(repodir,'results','group_level',['within_group' num2str(i) '_robust.mat']));
write(out.t, 'fname', fullfile(repodir,'results','group_level','cm_pre_to_post.nii'))
    
i=2; load(fullfile(repodir,'results','group_level',['within_group' num2str(i) '_robust.mat']));
write(out.t, 'fname', fullfile(repodir,'results','group_level','oxy_pre_to_post.nii'))

i=3; load(fullfile(repodir,'results','group_level',['within_group' num2str(i) '_robust.mat']));
write(out.t, 'fname', fullfile(repodir,'results','group_level','faml_pre_to_post.nii'))


%% full montage for between group comparisons 
o1 = canlab_results_fmridisplay([], 'montage_type', 'full', 'overlay', underlay);
o1 = multi_threshold(cm_vs_oxy, 'o2', o1, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);
saveas(gcf, fullfile(figdir, 'cm_vs_oxy_full_001_005'), 'jpeg')
%%
o1 = canlab_results_fmridisplay([], 'montage_type', 'full', 'overlay', underlay);
o1 = multi_threshold(cm_vs_faml, 'o2', o1, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);
saveas(gcf, fullfile(figdir, 'cm_vs_faml_full_001_005.tiff'))
%%
o1 = canlab_results_fmridisplay([], 'montage_type', 'full', 'overlay', underlay);
o1 = multi_threshold(oxy_vs_faml, 'o2', o1, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);
saveas(gcf, fullfile(figdir, 'oxy_vs_faml_full_001_005.jpeg'))

%%  full montage forfor each within-group comparison
for i=1:3
    load(fullfile(repodir,'results','group_level',['within_group' num2str(i) '_robust.mat']));
    out.t = apply_mask(out.t, gm_mask); % apply GM mask

    o1 = canlab_results_fmridisplay([], 'montage_type', 'full', 'overlay', underlay);
    o1 = multi_threshold(out.t, 'o2', o1, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);
    saveas(gcf, fullfile(figdir, ['within_group' num2str(i) '_full_001_005.png']))
end



%% ---------- slices for selective display ------------
% OK to use removeblobs() here b/c it works for slices


%% Fig 1a: Axial slices, group comparisons
% see also brainstem_slice3d. also surface_cutaway is amazing, but only shows medial surface

% set up montage
o1 = fmridisplay('overlay',underlay);
o1 = montage(o1, 'axial', 'slice_range', [-22 -16], 'onerow', 'spacing', 6, 'noverbose');

brighten(.2)
enlarge_axes(gcf, .8);
           
%% add blobs: CM vs Oxy
o1 = removeblobs(o1);
o1 = multi_threshold(cm_vs_oxy, 'o2', o1, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);
saveas(gcf, fullfile(figdir, ['cm_vs_oxy_axial.png'])); % print image to file

%% add blobs: CM vs Faml
o1 = removeblobs(o1);
o1 = multi_threshold(cm_vs_faml, 'o2', o1, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);

% also add in blob from randomize results
clusters_img = fmri_data(fullfile(repodir, 'data/for_randomise/masks/ofc_nac_vta/CMvFAML', 'clusters_img.nii'));
o1 = addblobs(o1, clusters_img);

%saveas(gcf, fullfile(figdir, ['cm_vs_faml_axial.png'])); % print image to file

%% Fig 1a: Sagittal slices: group comparisons
% set up montage
o2 = fmridisplay('overlay',underlay);
o2 = montage(o2, 'sagittal', 'slice_range', [-18 6], 'onerow', 'spacing', 6, 'noverbose');
brighten(.2)
enlarge_axes(gcf, .8);

%% add blobs: CM vs Oxy
o2 = removeblobs(o2);
o2 = multi_threshold(cm_vs_oxy, 'o2', o2, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);
saveas(gcf, fullfile(figdir, ['cm_vs_oxy_sag.png'])); % save image to file

%% add blobs: CM vs Faml
o2 = removeblobs(o2);
o2 = multi_threshold(cm_vs_faml, 'o2', o2, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);

% also add in blob from randomize results
clusters_img = fmri_data(fullfile(repodir, 'data/for_randomise/masks/ofc_nac_vta/CMvFAML', 'clusters_img.nii'));
o2 = addblobs(o2, clusters_img);

%saveas(gcf, fullfile(figdir, ['cm_vs_faml_sag.png'])); % save image to file




%% Fig 1b: Axial slices, effect within group
o3 = fmridisplay('overlay',underlay); brighten(.2);
o3 = montage(o3, 'axial', 'slice_range', [-16 4], 'onerow', 'spacing', 10, 'noverbose');
%enlarge_axes(gcf, 1.1) 

%% plot each group

for wh_image = 1:3 % 1 2 or 3, for group desired
    load(fullfile(repodir,'results','group_level',['within_group' num2str(wh_image) '_robust.mat']));
    out.t = apply_mask(out.t, gm_mask); % apply GM mask

    o3 = removeblobs(o3);
    o3 = multi_threshold(out.t, 'o2', o3, 'thresh', [.001 .005], 'sizethresh', [10 1]); %, 'poscolors', poscolor);

    % print image to file
    saveas(gcf, fullfile(figdir, ['within_group' num2str(wh_image) '_axial_001_005.png']))
end

%% which brainstem structure is it?  look at brainstem atlas Shen et al 2013

brainstem = load_atlas('brainstem') % Shen paper

i=[19 22 25];     % region 2 is closest
brainstem.labels{i}, brainstem.label_descriptions{i}
orthviews_multiple_objs({select_atlas_subset(brainstem, i), threshold(cm_vs_other, .005, 'unc')})

%% look at brainstem atlas Keuken

keuken = load_atlas('keuken') 
i=[1 3];   % left red nucleus, left SN
keuken.labels{i}, keuken.label_descriptions{i}
orthviews_multiple_objs({select_atlas_subset(keuken, i), threshold(cm_vs_other, .005, 'unc')})

%% or load in the PAG?
pag = canlab_load_ROI('pag') % i think hand-drawn from Tor
orthviews(pag)


%% table of results coordinates - cm vs control
dat_region = region(threshold(cm_vs_oxy, .001, 'unc', 'k', 10));
table(dat_region);

%% table of results coordinates - oxy vs faml
dat_region = region(threshold(cm_vs_faml, .001, 'unc', 'k', 10));
table(dat_region);

%% table of results coordinates - pre to post
for wh_image = 3 % 1 2 or 3, for group desired
    load(fullfile(repodir,'results','group_level',['within_group' num2str(wh_image) '_robust.mat']));
    out.t = apply_mask(out.t, gm_mask); % apply GM mask

    out.t = threshold(out.t, .001, 'unc', 'k', 10);
    dat_region = region(out.t);
    table(dat_region);
end

%% which anatomical regions are those Faml decreases in?
citi = load_atlas('cit168');
a = select_atlas_subset(citi, [16 5 8]);
a.label_descriptions
a.labels





%% interaction plots for ROIs that survive group x time

% load in full subj level data
time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')
grp = get_var(CM, 'Group');

%% CM vs Faml
rois = threshold(cm_vs_faml, .001, 'unc', 'k', 10);
roi_vals_T1 = extract_roi_averages(time1, rois, 'contiguous_regions');
roi_vals_T2 = extract_roi_averages(time2, rois, 'contiguous_regions');

for i=1:length(roi_vals)
    create_figure(['ROI ' num2str(i)]);
    
    h1 = lineplot_columns([roi_vals_T1(i).dat(grp==1) roi_vals_T2(i).dat(grp==1)], 'shade', 'color', 'r', 'x', [1 2], 'within');
    h1 = lineplot_columns([roi_vals_T1(i).dat(grp==3) roi_vals_T2(i).dat(grp==3)], 'shade', 'color', 'k', 'x', [1 2], 'within');
    h1 = lineplot_columns([roi_vals_T1(i).dat(grp==2) roi_vals_T2(i).dat(grp==2)], 'shade', 'color', 'b', 'x', [1 2], 'within');

    set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre-intervention' 'Post-intervention'}, 'FontSize', 24)
    title(['ROI center: ' num2str(roi_vals_T1(i).mm_center)])

end

%% examine buffering hypothesis

% load in all regions that showed a decrease in Faml
i=3; load(fullfile(repodir,'results','group_level',['within_group' num2str(i) '_robust.mat']));

% binarize
rois = threshold(out.t, .001, 'unc', 'k', 10);
%%
rois.dat = rois.sig;

% extract vals for all Ps in those regions
roi_vals_T1 = extract_roi_averages(time1, rois);
roi_vals_T2 = extract_roi_averages(time2, rois);

% make an interaction plot across the groups. don't sig test, just look at
% it descriptively, b/c that would be double dipping
create_figure('Faml decreases');
    
h1 = lineplot_columns([roi_vals_T1.dat(grp==1) roi_vals_T2.dat(grp==1)], 'shade', 'color', 'r', 'x', [1 2], 'within');
h2 = lineplot_columns([roi_vals_T1.dat(grp==3) roi_vals_T2.dat(grp==3)], 'shade', 'color', 'k', 'x', [1 2], 'within');
h3 = lineplot_columns([roi_vals_T1.dat(grp==2) roi_vals_T2.dat(grp==2)], 'shade', 'color', 'b', 'x', [1 2], 'within');

set(gca, 'XTick', [1 2], 'XTickLabel', {'Pre-intervention' 'Post-intervention'}, 'FontSize', 24)
legend([h1.line_han h3.line_han h2.line_han],{'CM' 'OxyPla' 'Faml'}, 'Location', 'SW')
saveas(gcf, fullfile(figdir, 'faml_decreases_roi_changes.svg'))