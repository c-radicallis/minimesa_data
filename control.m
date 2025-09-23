% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data2 in MATLAB with the System Identification Toolbox.
clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[6 60]}; opts1.PhaseMatching='on'; opts1.Grid='on';

% %% Data 1
% input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
% addpath(input_file_folder);
% file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
% file = 'pink_noise_40Hz_T3mm_scl=1_0.acq'; % load output acq
% scale = 1;
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);% load output acq
% n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
% data1 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data1.InputName  = 'x_drv_T_0';data1.OutputName = 'x_acq_T';data1.TimeUnit   = 'seconds';

%% data4
folder_18 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\1-8-2025\';
folder_58 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\5-8-2025\';
file = 'Noise_convertable_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder',folder_18 , 'OutputFolder', folder_18); % load input drv
file = 'noiseLTF_PID_10_0.1_0.1_run2_0.acq';
scale = 1;
LTF_to_TXT_then_load( file , 'InputFolder', folder_58 , 'OutputFolder', folder_58);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data4 = iddata(x_acq_T(1:nmin), scale*x_drv_T_0(1:nmin), Ts);data4.InputName  = 'x_drv_T_0';data4.OutputName = 'x_acq_T';data4.TimeUnit   = 'seconds';

%% Model training
nx = 3;
sys4 = n4sid(data4,nx,'Ts',Ts);% sys4.InputName  = data2.InputName; sys4.OutputName = data2.OutputName;

%% 
Ts_fpga= 1/1600;
z = tf('z',Ts_fpga);
Kc=10;
Ti=0.1;
Ki=Kc*Ts_fpga/Ti;
Controller = Kc + Ki/(z-1)
G_closed = d2d(sys4 , Ts_fpga)
% G_open = G_closed/(Controller*(1-G_closed)) 
G_open = minreal( G_closed/(Controller*(1-G_closed))) 
G_open.InputName =  'i_sv';

%%
tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
cutoff_frequency = 10; % Hz
PIDF  = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts)
G_PIDF_10Hz = feedback(PIDF*G_open, 1);

cutoff_frequency = 15; % Hz
PIDF   = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts);
G_PIDF_15Hz = feedback(PIDF*G_open, 1);

%%
obs = vpa(obsv(G_open));
r_obsv = rank(obs)
ctrlb = vpa(ctrb(G_open));
r_ctrlb = rank(ctrlb)

%
n_states=size(G_open.A,1);
G_open.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);

plant_aug = ss(G_open.A, G_open.B,[eye(n_states);G_open.C],G_open.D , Ts_fpga);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [G_open.StateName ; {'y_xT'}];  % plant output

sumblk1 = sumblk('e = x_ref - y_xT'); % Compute the error signal: e = r - y

integrator = tf(1,[1 0], Ts_fpga); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = diag([ones(1,n_states),1e2]);%1e3*diag([zeros(1,n_states),1]);%blkdiag(eye(nx), eye(ny));
R = eye(size(G_open.B,2));
K_lqi = lqi(G_open, Q, R)% Design the LQI controller for the original system

K  = K_lqi(1:n_states);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = [ G_open.StateName ; {'xi'}];
controller.OutputName = {'i_sv'};
Optimal_closed_loop = minreal(connect(plant_aug,  controller , integrator, sumblk1, 'x_ref','y_xT'));

figure;pzmap(Optimal_closed_loop);

%%
fig9 = figure(9);ax9 = axes(fig9); hold(ax9, 'on');
% bodeplot(sys4,opts1);
bodeplot(G_closed,opts1);
bodeplot(G_PIDF_10Hz,opts1);
bodeplot(G_PIDF_15Hz,opts1); 
bodeplot(Optimal_closed_loop,opts1);
legend();

%%
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
dados = load('elcentro.txt');
t_vector = dados(:,1);
ddx_ref = dados(:,2);
lim_displacement = 0.005; % m % Limits
lim_velocity = 0.4; % m/s
lim_force = 200e3; % N

s=tf('s') ;
scale=0.0375;
ddx_ref=scale*ddx_ref;
x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
% max_xref = max(x_ref);
% while max_xref > lim_displacement % Scaling down if necessary
%     scale = 0.95*scale;
%     ddx_ref = 0.95*ddx_ref;
%     x_ref = lsim(1/s^2,  ddx_ref , t_vector ,'foh');
%     max_xref = max(x_ref);
% end

%% Finding Response Spectre of Ground
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 2e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

[picos_ddx_ground , picos_x_ground] = ResponseSpectrum( t_vector , x_ref , ddx_ref, f_vector , 1);
x_cl = lsim(d2d(G_closed,Ts) ,  x_ref , t_vector,'zoh');
ddx_cl = secondDerivativeTime(x_cl , Ts);
[picos_ddx_cl , picos_x_cl] = ResponseSpectrum( t_vector , x_cl , ddx_cl, f_vector , 1);

x_PIDF_10Hz = lsim(d2d(G_PIDF_10Hz,Ts) ,  x_ref , t_vector,'zoh');
ddx_PIDF_10Hz = secondDerivativeTime(x_PIDF_10Hz , Ts);
[picos_ddx_PIDF_10Hz , picos_x_PIDF_10Hz] = ResponseSpectrum( t_vector , x_PIDF_10Hz , ddx_PIDF_10Hz, f_vector , 1);

x_PIDF_15Hz = lsim(d2d(G_PIDF_15Hz,Ts) ,  x_ref , t_vector,'zoh');
ddx_PIDF_15Hz = secondDerivativeTime(x_PIDF_15Hz , Ts);
[picos_ddx_PIDF_15Hz , picos_x_PIDF_15Hz] = ResponseSpectrum( t_vector , x_PIDF_15Hz , ddx_PIDF_15Hz, f_vector , 1);

x_Optimal = lsim(d2d(Optimal_closed_loop,Ts) ,  x_ref , t_vector,'zoh');
ddx_Optimal = secondDerivativeTime(x_Optimal , Ts);
[picos_ddx_Optimal , picos_x_Optimal] = ResponseSpectrum( t_vector , x_Optimal , ddx_Optimal, f_vector , 1);

baseFolder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\figures\';   % Base folder where you want to create the timestamped subfolder
ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
    mkdir(timeDir)
end

fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);%ylim([0 ceil(max(picos_ddx_T_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_ground,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_cl,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Standard -  MSE= %.2e',      mean((picos_ddx_ground-picos_ddx_cl).^2 )));% - Normal
plot(f_vector, picos_ddx_PIDF_10Hz,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'PIDF 10Hz - MSE= %.2e',   mean((picos_ddx_ground-picos_ddx_PIDF_10Hz).^2 )));
plot(f_vector, picos_ddx_PIDF_15Hz ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'PIDF 15Hz - MSE= %.2e', mean((picos_ddx_ground-picos_ddx_PIDF_15Hz).^2 )));
plot(f_vector, picos_ddx_Optimal ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Optimal Control -  MSE= %.2e', mean((picos_ddx_ground-picos_ddx_Optimal).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_ground,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_cl,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Standard - MSE= %.2e',     mean((picos_x_ground-picos_x_cl).^2 )));%- Normal
plot(f_vector, picos_x_PIDF_10Hz,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'PIDF 10Hz - MSE= %.2e',  mean((picos_x_ground-picos_x_PIDF_10Hz).^2 )));
plot(f_vector, picos_x_PIDF_15Hz, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'PIDF 15Hz - MSE= %.2e', mean((picos_x_ground-picos_x_PIDF_15Hz).^2 )));
plot(f_vector, picos_x_Optimal, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Optimal Control - MSE= %.2e', mean((picos_x_ground-picos_x_Optimal).^2 )));

set(fig8, 'WindowState', 'maximized');
exportgraphics(fig8,fullfile(timeDir,'Response_Spectra_N.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');

%%
figure; hold on; grid on;
plot(t_vector, x_ref);
plot(t_vector, x_cl);
plot(t_vector, x_PIDF_15Hz);
plot(t_vector , x_Optimal)