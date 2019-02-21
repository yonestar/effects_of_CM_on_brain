%% set up and load data
repodir = '/Users/yoni/Repositories/effects_of_CM_on_brain';
cd(repodir);

figdir = fullfile(repodir, 'figures');

load data/CM_data_included_Ss_only.mat

colors = {[1 .2 .2] [.2 .2 1] [.4 .4 .4] [.7 .7 .7]};

wh = get_var(CM, 'deltadon') > -45 | get_var(CM, 'deltadon') < -46 ; % drop the P who misunderstood donation instructions

%% bar graph of 3 groups
plot_var(CM, 'deltadon', 'subjtype', 'Group', 'color', colors, '95CI', 'wh_keep', wh)%, 'noviolin', 'noind') %deltadon and deltaFAS
set(gcf, 'Position',[ 840  92   228   247])
%set(gca, 'FontName', 'Helvetica Neue')

% print image to file
saveas(gcf, fullfile(figdir, 'deltadon_bar.svg'));

plot_var(CM, 'deltaFAS', 'subjtype', 'Group', 'color', colors, '95CI', 'wh', wh)%, 'noviolin', 'noind') %deltadon and deltaFAS
set(gcf, 'Position',[ 840  92   228   247])

% print image to file
saveas(gcf, fullfile(figdir, 'deltaFAS_bar.svg'));



 
%% glm of group differences in pre-to-post change in don or CC, including 2 and 3 way interactions
%out = glm(D, 'wd1', { 'med_vs_other' 'oxy_vs_faml' 'SSESdiff'},D.wh_keep.to_keep);% { 'dayscompleted'},D.wh_keep.to_keep); %to control for time1, add 'ccd1_1' to cell array
out = glm(CM, 'deltaFAS', { 'med_vs_oxy' 'faml_vs_other'},wh);% { 'dayscompleted'},D.wh_keep.to_keep); %to control for time1, add 'ccd1_1' to cell array

tcrit = tinv(.975, out.stat.dfe);
%%
out.stats = out.stat;
for i=1:length(out.stats.t)
    fprintf(', t(%d) = %3.2f, p = %3.2f, 95%% CI [%3.2f, %3.2f]\n',out.stats.dfe, out.stats.t(i), out.stats.p(i), out.stats.beta(i) - tcrit * out.stats.se(i), out.stats.beta(i) + tcrit * out.stats.se(i));
end

%% 2sample ttests on group diffs
var = get_var(CM, 'Exp_pre', wh);
group = get_var(CM, 'Group',wh);
sample1 = var(group==2); 
sample2 = var(group==3); 
[~,p, ci, stats] = ttest2(sample1,sample2);
fprintf(', T(%d) = %3.2f, p = %3.3f, 95%% CI [%3.2f, %3.2f]\n',stats.df, stats.tstat, p, ci(1), ci(2))

%% cohen's D for difference between two samples
varname = 'ccd1';
sample1 = get_var(D, varname, D.wh_keep.allmed); 
sample2 = get_var(D, varname, D.wh_keep.oxy); 
sample3 = get_var(D, varname, D.wh_keep.faml);

[sample1, sample2] = padwithnan(sample1, sample2, 1);
[sample1, sample3] = padwithnan(sample1, sample3, 1);

fprintf('CM vs OxyPla d = %3.2f\n', cohensd(sample1,sample2));
fprintf('CM vs Familiarity d = %3.2f\n', cohensd(sample1,sample3));
fprintf('OxyPla vs Familiarity d = %3.2f\n', cohensd(sample2,sample3));
fprintf('CM vs combined d = %3.2f\n', cohensd(sample1, [sample2; sample3]));

%% mean and 95% CI for each condition 
varname = 'Exp_pre';
all = get_var(CM, varname); 
sample1 = get_var(CM, varname, CM.wh_keep.med); 
sample2 = get_var(CM, varname, CM.wh_keep.oxy); 
sample3 = get_var(CM, varname, CM.wh_keep.faml);

fprintf('\n\n')
[~,~,ci,stats] = ttest(all);
fprintf('M = %3.2f, 95%% CI [%3.2f, %3.2f]\n', nanmean(all), ci(1), ci(2))

fprintf('\n')
[~,~,ci,stats] = ttest(sample1);
fprintf('M = %3.2f, 95%% CI [%3.2f, %3.2f]\n', nanmean(sample1), ci(1), ci(2))

[~,~,ci,stats] = ttest(sample2);
fprintf('M = %3.2f, 95%% CI [%3.2f, %3.2f]\n', nanmean(sample2), ci(1), ci(2))

[~,~,ci,stats] = ttest(sample3);
fprintf('M = %3.2f, 95%% CI [%3.2f, %3.2f]\n', nanmean(sample3), ci(1), ci(2))

%% ES and CI within group
varname = 'deltadon';
sample1 = get_var(CM, varname, CM.wh_keep.med); 
sample2 = get_var(CM, varname, CM.wh_keep.oxy); 
sample3 = get_var(CM, varname, CM.wh_keep.faml);

es = mes(sample3, 0, 'hedgesg', 'nBoot', 10000);
fprintf('g = %3.2f, 95%% CI [%3.2f, %3.2f]\n', es.hedgesg, es.hedgesgCi(1), es.hedgesgCi(2))

%% ES and CI between groups, using mes toolbox + bootstrapping
varname = 'deltadon';
sample1 = get_var(CM, varname, CM.wh_keep.med); 
sample2 = get_var(CM, varname, CM.wh_keep.oxy); 
sample3 = get_var(CM, varname, CM.wh_keep.faml);

es = mes(sample1, sample2, 'hedgesg', 'nBoot', 100000);
fprintf('g = %3.2f, 95%% CI [%3.2f, %3.2f]\n', es.hedgesg, es.hedgesgCi(1), es.hedgesgCi(2))




%% exact CIs for cohen's D, following Odgaard et al. 2010 -- from Emotion 2016 paper
%compute t
[~,p,ci,stats]=ttest2(sample1,[sample2; sample3]);
ncptvals = nctinv([.025 .975], stats.df, stats.tstat);
d = ncptvals / sqrt( (length(sample1)*length(sample2)) / (length(sample1)+length(sample2)) );
fprintf('[%3.2f, %3.2f]\n', d(1),d(2));

%% anova on various outcomes

% expectations: Exp_pre
% 'dayscompleted'  'daysPerCorr'
[P,ANOVATAB,STATS] = anova1(get_var(CM, 'Exp_pre'), get_var(CM, 'Group'))