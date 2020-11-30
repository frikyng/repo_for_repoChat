clear all;
close all;
load('test_set_1.mat')
[all_fits, all_pks] = detect_and_fit_events(all_concat, global_timescale, true);