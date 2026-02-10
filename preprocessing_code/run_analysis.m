close all;
clear all;

Paths.Data     = '\\\\your_path_to_analysis_code\RawData\';
Paths.plots     = '\\your_path_to_analysis_code\result_plot\';

addpath ''\\your_path_to_analysis_code\''

Const = storeConst();

%% BRAIN task (keyboard tapping)

behav_analysis_BRAIN_Test_reordered(Const, Paths);  
merge_behave_with_demographic(Const, Paths);