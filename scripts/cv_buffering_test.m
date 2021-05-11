% Test buffering hypothesis with unbiased voxel selection procedure. Does
% one CV iteration. In a CV loop, find 5% most decreased voxels in Faml
% group, extract data for those voxels in the held-out data for all 3
% groups.
%
% Input: 
%    - fmri_data object with a set of images representing pre-to-post change for a number N of participants (e.g., delta_relisten)
%    - assignments of those N images to group. This can be the true labels
%    or the shuffled labels.
%    - k: number of folds (e.g., 5)
%
% Output:
%    - CV predictions, for all Ss / all groups
%    - map of 5% most decreased voxels in Faml (sum across the folds, with each
%    unit representing a "hit" in one fold)

function [test_cm, test_oxy, test_faml, decreases] = cv_buffering_test(delta_relisten, grp, kfold)

    % get the groups
    delta_cm = get_wh_image(delta_relisten, grp==1);
    delta_oxy = get_wh_image(delta_relisten, grp==2);
    delta_faml = get_wh_image(delta_relisten, grp==3);

    % create train and test sets
    cvinds_faml = cvpartition(18, 'KFold', kfold);
    cvinds_cm = cvpartition(21, 'KFold', kfold);
    cvinds_oxy = cvpartition(18, 'KFold', kfold);
    
    decreases = [];
    
    for j = 1:kfold

        % find 5% most decreased voxels in Faml train set
        delta_faml_train = get_wh_image(delta_faml, training(cvinds_faml,j));
        t = ttest(delta_faml_train);
        decreases5prc = threshold(t, [prctile(t.dat, 5) Inf], 'raw-outside');
        decreases5prc.dat = decreases5prc.sig; % make into a "mask"
        %decreases5prc = threshold(t, .01, 'unc');    %figure; montage(decreases5prc)

        % save the area with decreases from the first fold
        if isempty(decreases)
            decreases = decreases5prc;
        else
            %decreases = image_math(decreases, decreases5prc, 'cat');
        end
        
        %% extract ROI vals from held out Ss
        
        delta_faml_test = get_wh_image(delta_faml, test(cvinds_faml,j));
        delta_cm_test = get_wh_image(delta_cm, test(cvinds_cm,j));
        delta_oxy_test = get_wh_image(delta_oxy, test(cvinds_oxy,j));
        
        tmp = extract_roi_averages(delta_faml_test, decreases5prc);
        test_faml(test(cvinds_faml,j)) = tmp.dat;

        tmp = extract_roi_averages(delta_cm_test, decreases5prc);
        test_cm(test(cvinds_cm,j)) = tmp.dat;

        tmp = extract_roi_averages(delta_oxy_test, decreases5prc);
        test_oxy(test(cvinds_oxy,j)) = tmp.dat;
        

    end % end cross val

end % end repeated cross val
