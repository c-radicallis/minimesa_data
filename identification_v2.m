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
%scale = 1;
LTF_to_TXT_then_load_wSV( file , folder_2910 , 'OutputFolder', folder_2910);

x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

% [c_drv,lags_drv] = xcorr(x_drv_T_0,x_acq_T);
% [c_sv,lags_sv] = xcorr(sv2_acq,x_acq_T,'normalized');
% figure(91); stem(lags_drv,c_drv); figure(92); stem(lags_sv,c_sv);

data11_openloop = iddata(x_acq_T, sv2_acq, Ts);data11_openloop.InputName  = 'sv2_acq';data11_openloop.OutputName = 'x_acq_T';data11_openloop.TimeUnit   = 'seconds';

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data11_closedloop =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);data11_closedloop.InputName  = 'x_drv_T_0';data11_closedloop.OutputName = 'x_acq_T';data11_closedloop.TimeUnit   = 'seconds';

%% Open Loop
g_data11_openloop = spa(data11_openloop, 1000);

% Model training
nx = 6 ;
n4sid_data11_openloop = n4sid(data11_openloop,nx,'Ts',Ts); n4sid_data11_openloop.InputName  = data11_openloop.InputName;n4sid_data11_openloop.OutputName = data11_openloop.OutputName;

% Resampled to 1600 Hz
Ts_fpga= 1/1600;
fpga_n4sid_data11_openloop = d2d(n4sid_data11_openloop, Ts_fpga);

% Figures Open Loop
fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
h = bodeplot(g_data11_openloop   ,opts1, "*");
bodeplot(n4sid_data11_openloop ,opts1, "m--")
bodeplot(fpga_n4sid_data11_openloop  ,opts1, "r--" )
showConfidence(h,3); legend(); grid on

    
%% PIDF tuning
% Ts_fpga= 1/1600;
% fpga_n4sid_data11_openloop = d2d(n4sid_data11_openloop, Ts_fpga);

tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
 
% cutoff_frequency = 5; % Hz
% PIDF  = pidtune(fpga_n4sid_data11_openloop,'PIDF',cutoff_frequency*2*pi,tuner_opts)
% G_PIDF_5Hz = feedback(PIDF*fpga_n4sid_data11_openloop, 1);

PIDF  = pidtune(fpga_n4sid_data11_openloop,'PIDF',tuner_opts)
G_PIDF_tracking = feedback(PIDF*fpga_n4sid_data11_openloop, 1);

G_PIDF_true_tune = feedback(6*fpga_n4sid_data11_openloop, 1);



%% Closed Loop
g_data11_closedloop = spa(data11_closedloop, 1000);

% Model training
nx = 20 ;
n4sid_data11_closedloop = n4sid(data11_closedloop,nx,'Ts',Ts); n4sid_data11_closedloop.InputName  = data11_closedloop.InputName;n4sid_data11_closedloop.OutputName = data11_closedloop.OutputName;

% Figures Closed Loop
fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop');
h = bodeplot(g_data11_closedloop   ,opts1, "*");
showConfidence(h,3);

bodeplot(n4sid_data11_closedloop ,opts1, "m-")
bodeplot(G_PIDF_tracking , opts1 )
bodeplot(G_PIDF_true_tune , opts1)
legend(); grid on

% figure(3), hold on;
% plot(time_drv_0,x_drv_T_0)
% plot(time_acq , x_acq_T)