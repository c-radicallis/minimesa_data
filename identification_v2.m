% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data2 in MATLAB with the System Identification Toolbox.
clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[0.7 40]};opts1.PhaseWrapping="on";%opts1.PhaseWrappingBranch=0;%opts1.Ylim={[-40 10]};

%% New data
% Data 11
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv

folder_2910 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\29-10-2025\';
file = 'pink_noise_40Hz_T3mm_scl1_P6_I0_D0_Fmax_0.acq'; % load output acq
scale = 1;
LTF_to_TXT_then_load_wSV( file , folder_2910 , 'OutputFolder', folder_2910);

% n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
% data1 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data1.InputName  = 'x_drv_T_0';data1.OutputName = 'x_acq_T';data1.TimeUnit   = 'seconds';

%% Old data
% Data 1
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv

file = 'pink_noise_40Hz_T3mm_scl=1_0.acq'; % load output acq
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);% load output acq

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data1 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data1.InputName  = 'x_drv_T_0';data1.OutputName = 'x_acq_T';data1.TimeUnit   = 'seconds';

% Data 2
file = 'pink_noise_40Hz_T3mm_scl=1.2_0.acq'; % load output acq
scale = 1.2;
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data2 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data2.InputName  = data1.InputName;data2.OutputName = data1.OutputName;data2.TimeUnit   = data1.TimeUnit;

% Data 3 
folder_18 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\1-8-2025\';
file = 'Noise1to200Hz_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

folder_58 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\5-8-2025\';
file = 'noise1to200hzLTF_PID_10_0.1_0.1_0.acq'; % load output acq
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data3 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data3.InputName  = data1.InputName;data3.OutputName = data1.OutputName;data3.TimeUnit   = data1.TimeUnit;

% Data 4
file = 'Noise_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

file = 'noiseLTF_PID_10_0.1_0.1_run2_0.acq';
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data4 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = data1.InputName;data4.OutputName = data1.OutputName;data4.TimeUnit   = data1.TimeUnit;

%
data_all = [data1; data2; data3;data4];

% Model training
nx = 5 ;

g_data1 = spa(data1, 1000);
n4sid_sys1 = n4sid(data1,nx,'Ts',Ts); n4sid_sys1.InputName  = data2.InputName;n4sid_sys1.OutputName = data2.OutputName;
g_data2 = spa(data2, 1000);
n4sid_sys2 = n4sid(data2,nx,'Ts',Ts); n4sid_sys2.InputName  = data2.InputName; n4sid_sys2.OutputName = data2.OutputName;
g_data3 = spa(data3, 1000);
n4sid_sys3 = n4sid(data3,nx,'Ts',Ts); n4sid_sys3.InputName  = data2.InputName;n4sid_sys3.OutputName = data2.OutputName;
g_data4 = spa(data4, 1000);
n4sid_sys4 = n4sid(data4,nx,'Ts',Ts); n4sid_sys4.InputName  = data2.InputName; n4sid_sys4.OutputName = data2.OutputName;
g_data_all = spa(data_all , 1000);
n4sid_sys_all = n4sid(data_all,nx,'Ts',Ts); n4sid_sys_all.InputName  = data2.InputName; n4sid_sys_all.OutputName = data2.OutputName;

%% Figures
fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on');
%bodeplot(g_data0   ,opts1, "*");
% bodeplot(g_data1   ,opts1, "r--");
% bodeplot(g_data2   ,opts1, "g--")
% bodeplot(g_data3   ,opts1, "b--")
% bodeplot(g_data4   ,opts1, "c--")
h = bodeplot(g_data_all,opts1, "m--");

showConfidence(h,3);

% bodeplot(n4sid_sys1   ,opts1, "r-")
% bodeplot(n4sid_sys2   ,opts1, "g-")
% bodeplot(n4sid_sys3   ,opts1, "b-")
% bodeplot(n4sid_sys4   ,opts1, "c-")
bodeplot(n4sid_sys_all,opts1, "m-")

legend(); grid on

% figure(2); compare(data1,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) 
% figure(3); compare(data2,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) 
% figure(4); compare(data3,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) 
% figure(5); compare(data4,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) 




