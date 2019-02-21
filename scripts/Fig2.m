clear all

%%
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
figdir = fullfile(repodir, 'figures');
outdir = fullfile(repodir, 'results', 'group_level');
cd(repodir);

gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));

% load between group contrasts
load(fullfile(repodir,'results','group_level','cm_vs_oxy_faml_robust.mat'));
cm_vs_oxy = get_wh_image(out2.t, 1);
cm_vs_oxy.dat = cm_vs_oxy.dat * -1; % flip so it is cm vs oxy, instead of oxy vs cm
cm_vs_oxy = apply_mask(cm_vs_oxy, gm_mask); % apply GM mask

cm_vs_faml = get_wh_image(out2.t, 2);
cm_vs_faml.dat = cm_vs_faml.dat * -1; % flip so it is cm vs faml, instead of faml vs cm
cm_vs_faml = apply_mask(cm_vs_faml, gm_mask); % apply GM mask

time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));
delta_relisten = image_math(time2, time1, 'minus'); 

load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')

grp = get_var(CM, 'Group');

%% get stats / error bars for group comparisons across the 7 networks / subcortical

[hh, buckner_med, labels] = wedge_plot_by_atlas(get_wh_image(delta_relisten, grp==1), 'atlases', {'subcortical_rl'})
%pause(5);
[hh, buckner_oxy, labels] = wedge_plot_by_atlas(get_wh_image(delta_relisten, grp==2), 'atlases', {'subcortical_rl'})
%pause(5)
[hh, buckner_famly, labels] = wedge_plot_by_atlas(get_wh_image(delta_relisten, grp==3), 'atlases', {'subcortical_rl'}, 'montage')

%% compute avg of SN and VTA
inds = [7 9 11]
labels{1}(inds)

buckner_med{1}(:,end+1) = mean(buckner_med{1}(:,inds),2);
buckner_oxy{1}(:,end+1) = mean(buckner_oxy{1}(:,inds),2);
buckner_famly{1}(:,end+1) = mean(buckner_famly{1}(:,inds),2);


%% test all atlas regions individually
[h,p]=ttest2(buckner_med{1}, buckner_famly{1})

g = mes(buckner_med{1}, buckner_famly{1}, 'hedgesg', 'nBoot', 10000, 'confLevel',0.90);
g.hedgesgCi


%% do the test avg across all subcortical structures together

figure;
violinplot({mean(buckner_med{:},2), mean(buckner_famly{:},2)})

g=mes(mean(buckner_med{:},2), mean(buckner_famly{:},2), 'hedgesg', 'nBoot', 10000)

dm(:,1) = scale([ones(21,1); zeros(18,1)], 1);
%dm(:,2) = ones(39,1);
[b,stats]=robustfit(dm, [mean(buckner_med{:},2); mean(buckner_famly{:},2)])

%% buckner atlas wedge plot on results from robust reg.

% group comparisons. before running these lines, set line 254 of
% wedge_plot_by_atlas to myouterradius = .6517; 
        
[hh, buckner_cm_vs_faml] = wedge_plot_by_atlas(cm_vs_faml, 'atlases', {'buckner'})
saveas(gcf, fullfile(figdir, 'wedge_buckner_cm_vs_faml.svg'))

[hh, buckner_cm_vs_oxy] = wedge_plot_by_atlas(cm_vs_oxy, 'atlases', {'buckner'})
saveas(gcf, fullfile(figdir, 'wedge_buckner_cm_vs_oxy.svg'))

%% subcortical atlas wedge plot on results from robust reg.

% group comparisons. before running these lines, set line 254 of
% wedge_plot_by_atlas to myouterradius = 2.1832 --> max wedge size across
% comparisons

[~, subc_cm_vs_faml] = wedge_plot_by_atlas(cm_vs_faml, 'atlases', {'subcortical_rl'})
saveas(gcf, fullfile(figdir, 'wedge_subcortical_cm_vs_faml.svg'))

[~, subc_cm_vs_oxy] = wedge_plot_by_atlas(cm_vs_oxy, 'atlases', {'subcortical_rl'})
saveas(gcf, fullfile(figdir, 'wedge_subcortical_cm_vs_oxy.svg'))

%% plot the subcortical montage, manually save
hh = wedge_plot_by_atlas(cm_vs_oxy, 'atlases', {'subcortical_rl'}, 'montage')


%% within group
for wh_image = 1:3 % 1 2 or 3, for group desired
    load(fullfile(repodir,'results','group_level',['within_group' num2str(wh_image) '_robust.mat']));
    [hh, output_values_by_region, labels, atlas_obj, colorband_colors] = wedge_plot_by_atlas(out.t, 'atlases', {'buckner'})
    saveas(gcf, fullfile(figdir, ['wedge_within_group_' num2str(wh_image) '.svg']))
end

%% make brain slices for buckner atlas reference image -- could also send the 'montage' flag into wedge_plot_by_atlas

o1 = fmridisplay('overlay',underlay);
o1 = montage(o1, 'axial', 'wh_slice', [0 0 -18], 'onerow', 'noverbose');
enlarge_axes(gcf, .8);
montage(buckner, 'nosymmetric', 'o2', o1)
saveas(gcf, fullfile(figdir, 'buckner_axial.png'))
%%
o2 = fmridisplay('overlay',underlay);
o2 = montage(o2, 'sagittal', 'wh_slice', [6 0 0], 'onerow', 'noverbose');
enlarge_axes(gcf, .8);
montage(buckner, 'nosymmetric', 'o2', o2)
saveas(gcf, fullfile(figdir, 'buckner_sag.png'))
