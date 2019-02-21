%% set up and load data
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
cd(repodir);

load results/group_level/CM_vs_Faml_robust.mat
dat = get_wh_image(out.t,1);

% load subject level data
time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));

% time2 - time1
delta_relisten = image_math(time2, time1, 'minus'); 

% gray matter mask; included in CanlabCore and SPM's DARTEL toolbox
gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));


%% make the ROIs

dat = apply_mask(dat, gm_mask);
dat = threshold(dat, .001, 'unc', 'k', 10); % only get clusters with > 10 voxels
regions = region(dat)

cerebellum = regions(1);
parahip = regions(2);
mOFC = regions(3);  
bstem = regions(4);
orthviews(mOFC);

% write out mask to file
%sgACC = region2fmri_data(sgACC, dat);
%write(sgACC, 'fname', fullfile(repodir, 'results', 'ROIs', 'sgACCmask_CMvsFamlrobust_001unc.nii'));


%% extract ROI activity

sgACCvals = extract_roi_averages(delta_relisten, region2fmri_data(mOFC, dat))
parahipvals = extract_roi_averages(delta_relisten, region2fmri_data(parahip, dat))
bstemvals = extract_roi_averages(delta_relisten, region2fmri_data(bstem, dat))
cerebvals = extract_roi_averages(delta_relisten, region2fmri_data(cerebellum, dat))

avg_all_ROIs = sum([sgACCvals.dat parahipvals.dat bstemvals.dat],2);

% add to CMdataset variable for easy plotting
load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat')) % load behav data
CM = add_vars(CM, [sgACCvals.dat parahipvals.dat bstemvals.dat cerebvals.dat avg_all_ROIs], 'Subj_Level', 'names', {'sgACC_change', 'parahip_change', 'brainstem_change' 'cereb_change' 'allROIs_change'});

%% make plot for brain changes, by group
plot_var(CM, 'sgACC_change', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind', 'nostars')
ylabel(''), xlabel(''), set(gca, 'XTickLabel',[]); set(gcf, 'Position', [345   384   209   298]);
saveas(gcf, fullfile(repodir, 'results', 'figures', 'sgACC_change.png'));
%%
plot_var(CM, 'parahip_change', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind', 'nostars')
ylabel(''), xlabel(''), set(gca, 'XTickLabel',[]); set(gcf, 'Position', [345   384   209   298]);
saveas(gcf, fullfile(repodir, 'results', 'figures', 'parahip_change.png'));

plot_var(CM, 'brainstem_change', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind', 'nostars')
ylabel(''), xlabel(''), set(gca, 'XTickLabel',[]); set(gcf, 'Position', [345   384   209   298]);
saveas(gcf, fullfile(repodir, 'results', 'figures', 'brainstem_change.png'));

plot_var(CM, 'cereb_change', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind')
plot_var(CM, 'allROIs_change', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind')

% also plot behav results
plot_var(CM, 'deltaFAS', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind')
ylabel(''), xlabel(''), set(gca, 'XTickLabel',[]); set(gcf, 'Position', [345   384   209   298]);
saveas(gcf, fullfile(repodir, 'results', 'figures', 'FASchange.png'));

plot_var(CM, 'deltadon', 'subjtype', 'Group', 'color', colors, '95CI', 'noviolin', 'noind')
ylabel(''), xlabel(''), set(gca, 'XTickLabel',[]); set(gcf, 'Position', [345   384   209   298]);
saveas(gcf, fullfile(repodir, 'results', 'figures', 'Donchange.png'));


%% effect of deltaFAS on ROI changes controling for group

% set group intercepts
g = get_var(CM, 'Group');
x = zeros(55,2);
x(g==1, 1) = 1;
x(g==2, 2) = 1;

% remove group effects from deltaFAS and deltadon
fitglm([get_var(CM, 'deltaFAS') x], avg_all_ROIs, 'VarNames', {'deltaFAS' 'CM' 'Oxy' 'ROIchange'})

%% effect of deltadon on ROI changes controling for group

% when correlating w/ don, drop the one P who misunderstood donation instructions
deltadon = get_var(CM, 'deltadon');
deltadon(49) = [];
x2 = x;
x2(49,:) = [];
avg_all_ROIs2 = avg_all_ROIs;
avg_all_ROIs2(49) = [];

fitglm([deltadon x2], avg_all_ROIs2, 'VarNames', {'deltadon' 'CM' 'Oxy' 'ROIchange'})
