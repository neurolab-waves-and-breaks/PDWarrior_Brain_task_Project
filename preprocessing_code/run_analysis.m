close all;
clear all;

Paths.Data     = '\\rdsfcifs.acrc.bris.ac.uk\Fischer_StudentProjects\PhD_Rui_Ni\Brain_Task_RawData\';
Paths.plots     = 'D:\PD_warrior_data\PD_warrior_tests\result_plot\';

%Paths.plots
% Paths.matfiles = 'D:\PD_warrior_tests\matfiles\FORCE_task\';
% Paths.plots    = 'C:\Users\ti21392\OneDrive - University of Bristol\Bristol\Experiments\PD_warrior\plots\';

% REQUIREMENTS: Need to download and add to the path: http://download.ant-neuro.com/matlab/
addpath('D:\New Volume\analysis_code_Rui\analysis_code\libeep-3.3.177-matlab')
addpath 'D:\New Volume\analysis_code_Rui\analysis_code'
addpath 'D:\'
Const = storeConst();

%% BRAIN task (keyboard tapping)

behav_analysis_BRAIN_Test_reordered(Const, Paths)  
