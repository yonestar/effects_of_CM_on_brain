% help find_closest_cluster
% help cluster_graphic_select

%% load
glasser = load_atlas('glasser');
r = atlas2region(glasser);
whos r
orthviews(r, 'unique');

%% click to cluster

%% get ID
xyz = spm_orthviews('Pos')
[clusters,wh] = find_closest_cluster(r, xyz)

% then find label, and select_atlas_subset using that label below

%orthviews(clusters);
%help cluster_graphic_select
%[clout,cl] = cluster_graphic_select(r);


%% dmPFC
dmPFC = select_atlas_subset(glasser, {'9m'});
orthviews(dmPFC)

%% mPFC

mPFC = select_atlas_subset(glasser, {'Ctx_10'}, 'flatten');
orthviews(mPFC)
