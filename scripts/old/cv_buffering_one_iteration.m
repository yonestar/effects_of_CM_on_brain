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



%% create train and test sets

kfold = 5; % doesn't seem to make a big diff whether k = 18 or k = 5, tho did not test thoroughly

[test_cm, test_faml, test_oxy] = deal([]);

for i=1:100 % number of times to repeat the CV loop, for getting a better point estimate of true CV effect size

    i

    cvinds_faml = crossvalind('Kfold',18,kfold);
    cvinds_cm = crossvalind('Kfold',21,kfold);
    cvinds_oxy = crossvalind('Kfold',18,kfold);

    for j = 1:kfold

        j 

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

        decreases5prc.dat = decreases5prc.sig; % make into a "mask"

        tmp = extract_roi_averages(delta_faml_test, decreases5prc);
        test_faml(i,cvinds_faml==j) = tmp.dat;

        tmp = extract_roi_averages(delta_cm_test, decreases5prc);
        test_cm(i,cvinds_cm==j) = tmp.dat;

        tmp = extract_roi_averages(delta_oxy_test, decreases5prc);
        test_oxy(i,cvinds_oxy==j) = tmp.dat;


    end % end cross val

end % end repeated cross val

%% save output



%% estimate buffering ES from the CV loop
for i=1:100
    g = mes(test_cm(i,:), test_faml(i,:), 'hedgesg', 'nBoot', 1000); 
    gs(i) = g.hedgesg;
end


%% compare mean pre-to-post changes in the CV-selected ROI across groups, for each CV repetition



figure; histogram(mean(test_cm') - mean(test_faml'),5)


%%

%save(fullfile(repodir,'data','cv_buffering_5prc_splithalf.mat'), 'cm_vs_faml', 'oxy_vs_faml', 'cm_vs_oxy')

% plot distributions

figure
violinplot({test_cm test_oxy test_faml},10)

% test





%g2 = mes(oxy_vs_faml(:), 0, 'hedgesg', 'nBoot', 1000)

[~,p]=ttest2(test_cm, test_faml)