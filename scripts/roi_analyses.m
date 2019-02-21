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

colors = {[1 .2 .2] [.2 .2 1] [.4 .4 .4] [.7 .7 .7]};

%% load atlases
subcortical = load_atlas('subcortical_rl');
glasser = load_atlas('glasser');


%% make my atlas with regions I want

% BA 10 = mPFC
combined_atlas = merge_atlases(select_atlas_subset(glasser, {'Ctx_10'}, 'flatten'), select_atlas_subset(glasser, {'OFC'}, 'flatten'));
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(subcortical, {'NAC'}));
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(subcortical, {'VTA'}));

% dmPFC (9m) and TPJ (PFm
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(glasser, {'9m'}, 'flatten'));
combined_atlas = merge_atlases(combined_atlas, select_atlas_subset(glasser, {'PFm'}, 'flatten'));

orthviews(combined_atlas);

%% extract from all ROIs for the 3 groups

for i = 1:3 % for the 3 groups  (inefficient code -- faster to extract all together and then divide into groups, but whatevs
    
    % DA pathway: extract NAC, VTA, mPFC, OFC
    vals = extract_roi_averages(get_wh_image(delta_relisten, grp==i), combined_atlas);
    mPFC{i} = vals(1).dat;
    OFC{i} = vals(2).dat;
    nac{i} = vals(3).dat;
    vta{i} = vals(4).dat;
   
    dmPFC{i} = vals(5).dat;
    TPJ{i} = vals(6).dat;
end

%% plot and test
clc

% choose var
var = OFC; roi_name = 'mOFC';

create_figure(roi_name); barplot_columns(var, roi_name, [], 'nofig', 'color', colors, 'noviolin', 'noind')
%set(gca, 'XTickLabels', {'CM' 'OxyPla' 'Faml'}, 'YTickLabels',[],'FontSize', 18); xlabel([]), ylabel([]), title(roi_name)
set(gca, 'XTickLabels', [], 'YTickLabels',[],'FontSize', 18); xlabel([]), ylabel([]), title(roi_name)
set(gcf, 'Position',[ 731   251   200   215])
saveas(gcf, fullfile(figdir, ['roi_' roi_name '.svg']))

%% compute group diffs
var = TPJ; clc

g = mes(var{1}, var{2}, 'hedgesg', 'nboot', 10000, 'confLevel', .9); 
fprintf('g = %3.2f, 90%% CI [%3.2f %3.2f]\n', g.hedgesg, g.hedgesgCi(1), g.hedgesgCi(2)) % bootstrap

g = mes(var{1}, var{3}, 'hedgesg', 'nboot', 10000, 'confLevel', .9); 
fprintf('g = %3.2f, 90%% CI [%3.2f %3.2f]\n', g.hedgesg, g.hedgesgCi(1), g.hedgesgCi(2)) % bootstrap

%dm(:,1) = scale([ones(21,1); zeros(18,1)], 1); [b,stats]=robustfit(dm, [var{1}; var{3}]); stats.p(2) % robust GLM


%% plot reference brain: surface (not working w/ colors)

% try to get surface display -- but its not showing w/ different colors
figure;
o2 = fmridisplay;
o2 = surface(o2, 'direction', 'right', 'orn', 'lateral');


for i=1:1
    o2 = montage(select_atlas_subset(combined_atlas,i), 'o2', o2, 'MarkerFaceColor', 'b');
end

%% reference brain: sag slice

overlay=which('keuken_2014_enhanced_for_underlay.img');
o1 = fmridisplay('overlay',overlay);
o1 = montage(o1, 'saggital', 'wh_slice', [-4 0 0], 'noverbose');
o1 = montage(combined_atlas, 'o2', o1)
enlarge_axes(gcf, .8)

saveas(gcf, fullfile(figdir, 'roi_ref.jpg'))

%% reference brain: display regions in separate figure

 montage(combined_atlas, 'regioncenters');






%% OLD -- load Denny
dennydir = '/Users/yoni/Repositories/Neuroimaging_Pattern_Masks/CANlab_Meta_analysis_maps/2012_Denny_SOMA_self_other_meta';
self_other = fmri_data(fullfile(dennydir,'SelfOthe_Other-Self.img'));
%%
denny_regions = threshold(self_other, [-Inf .055], 'raw-outside', 'k', 200);

% smooth and binarize
denny_regions = preprocess(denny_regions, 'smooth', 3)
denny_regions.dat = denny_regions.dat > .02;

orthviews(denny_regions)
%%
denny_regions = region(denny_regions) % 1=rTPJ, 2=lTPJ, 3 = PCC, 4 = dmPFC
denny_regions(3) = []; % drop PCC, no strong prior hypotheses here