% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data2 in MATLAB with the System Identification Toolbox.
clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 40]};

%% Training data2
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
addpath(input_file_folder);
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv

file = 'pink_noise_40Hz_T3mm_scl=1_0.acq'; % load output acq
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);% load output acq

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data1 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data1.InputName  = 'x_drv_T_0';data1.OutputName = 'x_acq_T';data1.TimeUnit   = 'seconds';

%% Validation  data2
file = 'pink_noise_40Hz_T3mm_scl=1.2_0.acq'; % load output acq
scale = 1.2;
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data2 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data2.InputName  = data1.InputName;data2.OutputName = data1.OutputName;data2.TimeUnit   = data1.TimeUnit;

%% Validation  data2 2

folder_18 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\1-8-2025\';
file = 'Noise1to200Hz_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

folder_58 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\5-8-2025\';
file = 'noise1to200hzLTF_PID_10_0.1_0.1_0.acq'; % load output acq
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data3 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data3.InputName  = data1.InputName;data3.OutputName = data1.OutputName;data3.TimeUnit   = data1.TimeUnit;

%DOESNT SEEM TO BE RIGHT %  % file = 'noise1to200hzLTF_scale0.5_PID_10_0.1_0.1.acq'; % load output acq   
%DOESNT SEEM TO BE RIGHT %  file ='noiseLTF_scale0.5_PID_10_0.1_0.1_0.acq';
% MAY BE GENERALIZED ISSUE WITH SCALED DRIVERS FROM THE NOISE LTF's

%% data4
file = 'Noise_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

file = 'noiseLTF_PID_10_0.1_0.1_0.acq';
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data4 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = data1.InputName;data4.OutputName = data1.OutputName;data4.TimeUnit   = data1.TimeUnit;

%% val_data4

file = 'Noise_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv

file = 'noiseLTF_PID_10_0.1_0.1_run2_0.acq';
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);

n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data4 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = data1.InputName;data4.OutputName = data1.OutputName;data4.TimeUnit   = data1.TimeUnit;


%% Model training
nx = 3;

n4sid_sys1 = n4sid(data1,nx,'Ts',Ts);n4sid_sys1.InputName  = data2.InputName;n4sid_sys1.OutputName = data2.OutputName;
n4sid_sys2 = n4sid(data2,nx,'Ts',Ts); n4sid_sys2.InputName  = data2.InputName; n4sid_sys2.OutputName = data2.OutputName;
n4sid_sys3 = n4sid(data3,nx,'Ts',Ts); n4sid_sys3.InputName  = data2.InputName;n4sid_sys3.OutputName = data2.OutputName;
n4sid_sys4 = n4sid(data4,nx,'Ts',Ts); n4sid_sys4.InputName  = data2.InputName; n4sid_sys4.OutputName = data2.OutputName;

%% Figures
% fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on');
% bodeplot(n4sid_sys1,opts1)
% bodeplot(n4sid_sys2,opts1)
% bodeplot(n4sid_sys3,opts1)
% bodeplot(n4sid_sys4,opts1)
% legend(); grid on
% 
% figure(2); compare(data1,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) %title('Model Training');
% 
% figure(3); compare(data2,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) %title('Model Validation');
% 
% figure(4); compare(data3,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) %title('Model Validation');
% 
% figure(5); compare(data4,n4sid_sys1,n4sid_sys2,n4sid_sys3,n4sid_sys4) %title('Model Validation');

%% Proportional=1 &  I=D=0
% file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
% 
% file = 'pink_noise_40Hz_T3mm_0_P=1_I=D=0_scale0.5_0.acq';
% scale = 0.5;
% LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);
% 
% n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
% data5 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = data1.InputName;data4.OutputName = data1.OutputName;data4.TimeUnit   = data1.TimeUnit;
% 
% %
% file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
% file = 'pink_noise_40Hz_T3mm_0_P=1_I=D=0_0.acq';
% scale = 1;
% LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);
% 
% n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
% data6 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = data1.InputName;data4.OutputName = data1.OutputName;data4.TimeUnit   = data1.TimeUnit;
% 
% %
% file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
% file = 'pink_noise_40Hz_T3mm_0_P=1_I=D=0_scale1.5_0.acq';
% scale = 1.5;
% LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);
% 
% n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
% data7 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = data1.InputName;data4.OutputName = data1.OutputName;data4.TimeUnit   = data1.TimeUnit;
% 
% n4sid_sys5 = n4sid(data5,nx,'Ts',Ts); n4sid_sys5.InputName  = data2.InputName; n4sid_sys5.OutputName = data2.OutputName;
% n4sid_sys6 = n4sid(data6,nx,'Ts',Ts); n4sid_sys6.InputName  = data2.InputName; n4sid_sys6.OutputName = data2.OutputName;
% n4sid_sys7 = n4sid(data7,nx,'Ts',Ts); n4sid_sys7.InputName  = data2.InputName; n4sid_sys7.OutputName = data2.OutputName;
% 
% figure(6); compare(data5,n4sid_sys5,n4sid_sys6,n4sid_sys7) ; title('P=1 & I=D=0');
% figure(7); compare(data6,n4sid_sys5,n4sid_sys6,n4sid_sys7) ; title('P=1 & I=D=0');
% figure(8); compare(data7,n4sid_sys5,n4sid_sys6,n4sid_sys7) ; title('P=1 & I=D=0');
% 
% fig9 = figure(9);ax9 = axes(fig9); hold(ax9, 'on');
% bodeplot(n4sid_sys5,opts1)
% bodeplot(n4sid_sys6,opts1)
% bodeplot(n4sid_sys7,opts1)
% legend(); grid on

%% 
z = tf('z',Ts);
Kc=10;
Ti=0.1;
Ki=Kc*Ts/Ti;
Controller = Kc + Ki*(1+z^-1)
G_closed = n4sid_sys4;
G_open = minreal(G_closed/(Controller*(1-G_closed)))

tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
cutoff_frequency = 10; % Hz
PIDF   = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts);
G_closed_PIDF_10Hz = feedback(PIDF*G_open, 1);
cutoff_frequency = 15; % Hz
PIDF   = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts);
G_closed_PIDF_15Hz = feedback(PIDF*G_open, 1);

% obs = vpa(obsv(G_open));
% r_obsv = rank(obs)
% ctrlb = vpa(ctrb(G_open));
% r_ctrlb = rank(ctrlb)

n_states=size(G_open.A,1);
G_open.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);

plant_aug = ss(G_open.A, G_open.B,[eye(n_states);G_open.C],G_open.D);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [G_open.StateName ; {'x_acq_T'}];  % plant output
sumblk1 = sumblk('e = x_drv_T_0 - x_acq_T'); % Compute the error signal: e = r - y
integrator = tf(1,[1 0]); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = 1e19*diag([zeros(1,n_states),1]);%blkdiag(eye(nx), eye(ny));
R = 1e-99*eye(size(G_open.B,2));
K_lqi = lqi(G_open, Q, R)% Design the LQI controller for the original system

K  = K_lqi(1:n_states);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = [ G_open.StateName ; {'xi'}];
controller.OutputName = {'i_sv'};
Optimal_closed_loop = connect(plant_aug,  controller , integrator, sumblk1, 'x_drv_T_0','x_acq_T')

close all;
fig9 = figure(9);ax9 = axes(fig9); hold(ax9, 'on');
bodeplot(G_closed,opts1);
bodeplot(G_closed_PIDF_10Hz,opts1);
bodeplot(G_closed_PIDF_15Hz,opts1); 
bodeplot(Optimal_closed_loop,opts1);
legend(); grid on;

% figure(10); compare(data3,G_closed,G_closed_PIDF_10Hz,G_closed_PIDF_15Hz,Optimal_closed_loop) ;
% ylim([-1 1]*6e-3)