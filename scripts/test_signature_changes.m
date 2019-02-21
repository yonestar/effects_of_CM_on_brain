% before running this script, must run prep_4_apply_signatures_and_save.m
% this saves signature response values in the DAT object.

basedir = '/Users/yoni/Repositories/effects_of_CM_on_brain';

load(fullfile(basedir, 'results', 'image_names_and_setup.mat'));

distress = DAT.SIG_contrasts.raw.cosine_sim.Empathic_Dist;
care = DAT.SIG_contrasts.raw.cosine_sim.Empathic_Care;

group = DAT.BETWEENPERSON.group;

d{1} = distress.pre_to_post_change