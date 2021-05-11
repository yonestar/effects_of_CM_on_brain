%% set up and load data
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
cd(repodir);

% load results/group_level/CM_vs_Faml_robust.mat
% dat = get_wh_image(out.t,1);

% load subject level data
time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

% time2 - time1
delta_relisten = image_math(time2, time1, 'minus'); 

% load the ROI
ofc = fmri_data(fullfile(repodir, 'data', 'for_randomise', 'masks','ofc_nac_vta', 'CMvFAML','clusters_img.nii.gz'));

% behav data
load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat')) % load behav data

%% extract ROI activity

ofc_vals = extract_roi_averages(delta_relisten, ofc)

%% sanity check

figure; boxchart(get_var(CM, 'Group'), ofc_vals.dat)

%% test for corr within CM group

group = get_var(CM, 'Group');
deltafas = get_var(CM, 'deltaFAS');
deltadon = get_var(CM, 'deltadon');

% remove group effects from deltaFAS and deltadon

create_figure('scat FAS');
plot_correlation(ofc_vals.dat(group==1), deltafas(group==1), 'robust', 'xlabel', 'Δ OFC', 'ylabel', 'Δ FAS');
set(gcf,'Position',[   680   780   240   198]);
saveas(gcf, fullfile(repodir, 'figures', 'within_CM_scatter_FAS.png'))

create_figure('scat don');
plot_correlation(ofc_vals.dat(group==1), deltadon(group==1), 'robust', 'xlabel', 'Δ OFC', 'ylabel', 'Δ Donation');


%% interaction plot

ofc1 = extract_roi_averages(time1, ofc)
ofc2 = extract_roi_averages(time2, ofc)
create_figure('int')
h1=lineplot_columns( {ofc1.dat(group==1) ofc2.dat(group==1)}, 'color', 'r', 'x', [1 2]); hold on
h1=lineplot_columns( {ofc1.dat(group==2) ofc2.dat(group==2)}, 'color', 'b', 'x', [1 2]); hold on
h1=lineplot_columns( {ofc1.dat(group==3) ofc2.dat(group==3)}, 'color', 'k', 'x', [1 2]); hold on

%h2=lineplot_columns( {tabl.(outcome)(wh_T1 & wh_wl) tabl.(outcome)(wh_T2 & wh_wl)}, 'color', [.3 .3 .3], 'x', [1.05 2.05]);  hold on

%lineplot(repmat([1 2], 57, 1), [ofc1.dat ofc2.dat])
