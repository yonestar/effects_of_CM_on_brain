clear all

%% load
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
figdir = fullfile(repodir, 'figures');
outdir = fullfile(repodir, 'results', 'group_level');
cd(repodir);

gm_mask = fmri_data(which('dartel_spm12_mni152_gray_matter_mask.img'));

time1 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time1', '*', 'con_0003.img')));
time2 = fmri_data(filenames(fullfile('data', 'subject level', 'relisten_time2', '*', 'con_0007.img')));
delta_relisten = image_math(time2, time1, 'minus');
delta_relisten = apply_mask(delta_relisten, gm_mask);

% divide into the the three groups
load(fullfile(repodir, 'data', 'CM_data_included_Ss_only.mat'), 'CM')
grp = get_var(CM, 'Group');

delta_cm = get_wh_image(delta_relisten, grp==1);
delta_oxy = get_wh_image(delta_relisten, grp==2);
delta_faml = get_wh_image(delta_relisten, grp==3);


% repeated CV, to get a better point estimate of true effect

for i = 1:20
    
    fprintf('Iteration %d\n', i);
    
    %% create train and test sets (split half)
    kfold = 5;
    cvinds_faml = crossvalind('Kfold',18,kfold);
    cvinds_cm = crossvalind('Kfold',21,kfold);
    cvinds_oxy = crossvalind('Kfold',18,kfold);
    
    for j = 1:kfold
        
        delta_faml_train = get_wh_image(delta_faml, cvinds_faml~=j);
        delta_cm_train = get_wh_image(delta_cm, cvinds_cm~=j);
        delta_oxy_train = get_wh_image(delta_oxy, cvinds_oxy~=j);
        
        delta_faml_test = get_wh_image(delta_faml, cvinds_faml==j);
        delta_cm_test = get_wh_image(delta_cm, cvinds_cm==j);
        delta_oxy_test = get_wh_image(delta_oxy, cvinds_oxy==j);
        
        %% find 5% most decreased voxels in Faml train set
        t = ttest(delta_faml_train);
        decreases5prc = threshold(t, [prctile(t.dat, 5) Inf], 'raw-outside');
        %decreases5prc = threshold(t, .01, 'unc');
        %figure; montage(decreases5prc)
        
        %% extract ROI vals from held out Ss
        decreases5prc.dat = decreases5prc.sig;
        test_faml = extract_roi_averages(delta_faml_test, decreases5prc);
        test_cm = extract_roi_averages(delta_cm_test, decreases5prc);
        test_oxy = extract_roi_averages(delta_oxy_test, decreases5prc);
        
        %% test group diffs
        
        %figure; violinplot({test_cm.dat test_oxy.dat test_faml.dat})
        
        % find mean diffs -- OR SHOULD I COMPUTE THE EFFECT SIZE HERE, AND
        % SAVE THAT?
        cm_vs_faml(i,j) = mean(test_cm.dat) - mean(test_faml.dat);
        oxy_vs_faml(i,j) = mean(test_oxy.dat) - mean(test_faml.dat);
        cm_vs_oxy(i,j) = mean(test_cm.dat) - mean(test_oxy.dat);
        
    end % end split half cross val
    
end % end iteration


%% save output

save(fullfile(repodir,'data','cv_buffering_5prc_5fold.mat'), 'cm_vs_faml', 'oxy_vs_faml', 'cm_vs_oxy')
%% plot distributions

figure
histogram(cm_vs_faml,10)

figure
histogram(oxy_vs_faml,10)

mean(cm_vs_faml)
mean(oxy_vs_faml)