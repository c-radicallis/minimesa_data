% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output val_data in MATLAB with the System Identification Toolbox.
clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);

%% Training val_data
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
addpath(input_file_folder);
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv

file = 'pink_noise_40Hz_T3mm_scl=1_0.acq'; % load output acq
scale = 1;

LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);% load output acq
Ts = time_acq(2);
train_data = iddata(x_acq_T, scale*x_drv_T_0, Ts);
train_data.InputName  = 'x_drv_T_0';
train_data.OutputName = 'x_acq_T';
train_data.TimeUnit   = 'seconds';

%% Validation  val_data
file = 'pink_noise_40Hz_T3mm_scl=1.2_0.acq'; % load output acq
scale = 1.2;

LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);

% n1 = numel(x_drv_T_0);
% n2 = numel(x_acq_T);
% nmin = min(n1, n2);
% x_drv_T_0 = x_drv_T_0(1:nmin); % Truncate both to the shortest length
% x_acq_T = x_acq_T(1:nmin);
% ddx_acq_T = ddx_acq_T(1:nmin);

val_data = iddata(x_acq_T, scale*x_drv_T_0, Ts);
val_data.InputName  = train_data.InputName;
val_data.OutputName = train_data.OutputName;
val_data.TimeUnit   = train_data.TimeUnit;

%% Validation  val_data 2

folder_18 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\1-8-2025\';
file = 'Noise1to200Hz_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

folder_58 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\5-8-2025\';
file = 'noise1to200hzLTF_PID_10_0.1_0.1_0.acq'; % load output acq
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);
n2 = numel(x_acq_T);
nmin = min(n1, n2);
x_drv_T_0 = x_drv_T_0(1:nmin); % Truncate boh to the shortest length
x_acq_T = x_acq_T(1:nmin);

val_data2 = iddata(x_acq_T, scale*x_drv_T_0, Ts);
val_data2.InputName  = train_data.InputName;
val_data2.OutputName = train_data.OutputName;
val_data2.TimeUnit   = train_data.TimeUnit;

%DOESNT SEEM TO BE RIGHT %  % file = 'noise1to200hzLTF_scale0.5_PID_10_0.1_0.1.acq'; % load output acq   

% val_data3

file = 'Noise_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

%DOESNT SEEM TO BE RIGHT %  file ='noiseLTF_scale0.5_PID_10_0.1_0.1_0.acq';
% MAY BE GENERALIZED ISSUE WITH SCALED DRIVERS FROM THE NOISE LTF's
scale = 0.5;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);
n2 = numel(x_acq_T);
nmin = min(n1, n2);
x_drv_T_0 = x_drv_T_0(1:nmin); % Truncate boh to the shortest length
x_acq_T = x_acq_T(1:nmin);

val_data3 = iddata(x_acq_T, scale*x_drv_T_0, Ts);
val_data3.InputName  = train_data.InputName;
val_data3.OutputName = train_data.OutputName;
val_data3.TimeUnit   = train_data.TimeUnit;

%% Model training
nx = 3;

n4sid_sys = n4sid(train_data,nx,'Ts',Ts)
n4sid_sys.InputName  = val_data.InputName;
n4sid_sys.OutputName = val_data.OutputName;

n4sid_sys_1 = n4sid(val_data,nx,'Ts',Ts)
n4sid_sys_1.InputName  = val_data.InputName;
n4sid_sys_1.OutputName = val_data.OutputName;

n4sid_sys_2 = n4sid(val_data2,nx,'Ts',Ts)
n4sid_sys_2.InputName  = val_data.InputName;
n4sid_sys_2.OutputName = val_data.OutputName;

n4sid_sys_3 = n4sid(val_data3,nx,'Ts',Ts)
n4sid_sys_3.InputName  = val_data.InputName;
n4sid_sys_3.OutputName = val_data.OutputName;

%% Figures
fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on');
opts1=bodeoptions('cstprefs');
opts1.FreqUnits = 'Hz';
opts1.XLim={[1 40]};
bodeplot(n4sid_sys,opts1)
bodeplot(n4sid_sys_1,opts1)
bodeplot(n4sid_sys_2,opts1)
bodeplot(n4sid_sys_3,opts1)
legend()
grid on

figure(2); 
compare(train_data,n4sid_sys,n4sid_sys_1,n4sid_sys_2,n4sid_sys_3)
%title('Model Training');

figure(3); 
compare(val_data,n4sid_sys,n4sid_sys_1,n4sid_sys_2,n4sid_sys_3) 
%title('Model Validation');

figure(4); 
compare(val_data2,n4sid_sys,n4sid_sys_1,n4sid_sys_2,n4sid_sys_3) 
%title('Model Validation');

figure(5); 
compare(val_data3,n4sid_sys,n4sid_sys_1,n4sid_sys_2,n4sid_sys_3) 
%title('Model Validation');