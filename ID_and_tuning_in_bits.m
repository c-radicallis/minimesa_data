clear;clc;close all; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-left
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
Ts_fpga= 1/5000;

% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;%opts1.Ylim={[-40 10]};

% Control channel AI2 Displacement - 16 bit signed integer to mm conversion
a = 0.000485;
b = -0.2;
bits2mm = @(bits) a*bits+b;
mm2bits = @(mm) (mm-b)/a;
clear a b

fir_np=100;
np_CL=4;
np_OL=4;

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = mm2bits(x_drv_T_0*1e3); % convert to mm
clear x_drv_L_0  x_drv_V_0

% %%  data_P7
% folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
% file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
% x_acq_T = mm2bits(x_acq_T*1e3);
% sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped
% Kp=7
% results_P7_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);
% 
% %%  Data P10
% folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
% file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
% x_acq_T = mm2bits(x_acq_T*1e3);
% sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped
% Kp=10
% results_P10_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = mm2bits(x_acq_T*1e3);
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped
Kp=15
results_P15_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data2 in MATLAB with the System Identification Toolbox.
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[5 100]}; opts1.PhaseMatching='on'; opts1.Grid='on';

%clc;
%close all;

OL_direct = results_P15_pink.OL_direct;
OL_indirect = results_P15_pink.OL_indirect;
OL_est_nonLin = ss(results_P15_pink.OL_est_nonLin);
CL_from_OL_direct = results_P15_pink.CL_from_OL_direct;
CL = results_P15_pink.CL;
CL_from_OL_est_nonLin = results_P15_pink.CL_from_OL_est_nonLin;

% n=0;
% while n<2
% figure(1);hold on;
% pzmap(OL_indirect)
% legend;
% figure(2);hold on;
% bodeplot(OL_indirect,opts1)
% legend;
% OL_indirect =minreal(results_P15_pink.OL_indirect)
% n=n+1;
% end
% 
% figure;hold on;
% pzmap(OL_direct,OL_indirect,OL_est_nonLin)
% legend;

OL_indirect =minreal(results_P15_pink.OL_indirect)

%
G_open=d2d(OL_est_nonLin,Ts_fpga,'tustin');
G_closed = d2d(results_P15_pink.CL_from_OL_est_nonLin,Ts_fpga,'tustin');

% figure;hold on;bodeplot(OL_est_nonLin,G_open,opts1);legend;

tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');

cutoff_frequency = 10; % Hz
PIDF  = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts)
G_PIDF_10Hz = feedback(PIDF*G_open, 1);
cutoff_frequency = 15; % Hz
PIDF   = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts)
G_PIDF_15Hz = feedback(PIDF*G_open, 1);
cutoff_frequency = 20; % Hz
PIDF  = pidtune(G_open,'PIDF',cutoff_frequency*2*pi,tuner_opts)
G_PIDF_20Hz = feedback(PIDF*G_open, 1);

%
% obs = vpa(obsv(G_open));
% r_obsv = rank(obs)
% ctrlb = vpa(ctrb(G_open));
% r_ctrlb = rank(ctrlb)

n_states=size(G_open.A,1);
G_open.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);

plant_aug = ss(G_open.A, G_open.B,[eye(n_states);G_open.C],0 , Ts_fpga);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [G_open.StateName ; {'y_xT'}];  % plant output

sumblk1 = sumblk('e = x_ref - y_xT'); % Compute the error signal: e = r - y

integrator = tf(1,[1 0], Ts_fpga); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = diag([0*ones(1,n_states),1e2]);%1e3*diag([zeros(1,n_states),1]);%blkdiag(eye(nx), eye(ny));
R = eye(size(G_open.B,2));
K_lqi = lqi(G_open, Q, R)% Design the LQI controller for the original system

K  = K_lqi(1:n_states);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = [ G_open.StateName ; {'xi'}];
controller.OutputName = {'i_sv'};
Optimal_CL = connect(plant_aug,  controller , integrator, sumblk1, 'x_ref','y_xT');
% Optimal_CL_minreal = minreal(Optimal_CL,1e-4)
% 
% figure;hold on;pzmap(Optimal_CL,Optimal_CL_minreal);legend;

figure;hold on;
bodeplot(G_closed,opts1);
bodeplot(G_PIDF_10Hz,opts1);
bodeplot(G_PIDF_15Hz,opts1); 
bodeplot(G_PIDF_20Hz,opts1);    
bodeplot(Optimal_CL,opts1);
legend;

%%  Data wc_10Hz
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\pink_noise';
file = 'wc_10Hz_0.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201 , 'OutputFolder', folder_1201);
x_acq_T = mm2bits(x_acq_T*1e3);
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped
% w_c = 10Hz
Kp = 6.75;
Ki = 68;
Kd = .0712891;
Tf = 0.00161743;
controller_at_200 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
results_wc10_pink = twoStageMethod( controller_at_200 , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%%  Data wc_15Hz
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\pink_noise';
file = 'wc_15Hz_0.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201 , 'OutputFolder', folder_1201);
x_acq_T = mm2bits(x_acq_T*1e3);
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped
% w_c = 15Hz
Kp = 8.75
Ki = 139
Kd = 0.125977
Tf = 0.00109863
controller_at_200 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
results_wc15_pink = twoStageMethod(controller_at_200 , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);


%%  Data wc_20Hz
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\pink_noise';
file = 'wc_20Hz_0.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201 , 'OutputFolder', folder_1201);
x_acq_T = mm2bits(x_acq_T*1e3);
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped
% w_c = 20Hz
Kp = 8.562500;
Ki = 75;
Kd = 0.18457;
Tf = 0.000167847;
controller_at_200 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
results_wc20_pink = twoStageMethod(controller_at_200 , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);




%% Let's compare open loop tranfers fucntions
close all;
%opts1.MagLowerLim = -60; opts1.MagLowerLimMode='manual';
opts1.XLim={[0.1 50]};

figure; hold on; opts1.Title.String='OL_{direct}';
% bodeplot(results_P7_sineSweep.OL_direct, opts1);
% bodeplot(results_P15_sineSweep.OL_direct, opts1);
% bodeplot(results_P5_pink.OL_direct, opts1);
bodeplot(results_P7_pink.OL_direct, opts1);
bodeplot(results_P10_pink.OL_direct, opts1);
bodeplot(results_P15_pink.OL_direct, opts1);
bodeplot(results_wc10_pink.OL_direct, opts1);
bodeplot(results_wc15_pink.OL_direct, opts1);
bodeplot(results_wc20_pink.OL_direct, opts1);
legend('P7','P10','P15','wc10','wc15','wc20'); grid on;

figure; hold on; opts1.Title.String='OL_{indirect}';
% bodeplot(results_P7_sineSweep.OL_indirect, opts1);
% bodeplot(results_P15_sineSweep.OL_indirect, opts1);
% bodeplot(results_P5_pink.OL_indirect, opts1);
bodeplot(results_P7_pink.OL_indirect, opts1);
bodeplot(results_P10_pink.OL_indirect, opts1);
bodeplot(results_P15_pink.OL_indirect, opts1);
bodeplot(results_wc10_pink.OL_indirect, opts1);
bodeplot(results_wc15_pink.OL_indirect, opts1);
bodeplot(results_wc20_pink.OL_indirect, opts1);
legend('P7','P10','P15','wc10','wc15','wc20'); grid on;

figure; hold on; opts1.Title.String='OL_{est nonLin}';
% bodeplot(results_P7_sineSweep.OL_est_nonLin, opts1);
% bodeplot(results_P15_sineSweep.OL_est_nonLin, opts1);
% bodeplot(results_P5_pink.OL_est_nonLin, opts1);
bodeplot(results_P7_pink.OL_est_nonLin, opts1);
bodeplot(results_P10_pink.OL_est_nonLin, opts1);
bodeplot(results_P15_pink.OL_est_nonLin, opts1);
bodeplot(results_wc10_pink.OL_est_nonLin, opts1);
bodeplot(results_wc15_pink.OL_est_nonLin, opts1);
bodeplot(results_wc20_pink.OL_est_nonLin, opts1);
legend('P7','P10','P15','wc10','wc15','wc20'); grid on;

% Let's compare CLOSED LOOP tranfers fucntions
%opts1.MagLowerLim = -40;
opts1.XLim={[5 50]};

figure; hold on; opts1.Title.String='CL from OL_{direct}';
% bodeplot(results_P7_sineSweep.CL_from_OL_direct, opts1);
% bodeplot(results_P15_sineSweep.CL_from_OL_direct, opts1);
% bodeplot(results_P5_pink.CL_from_OL_direct, opts1);
bodeplot(results_P7_pink.CL_from_OL_direct, opts1);
bodeplot(results_P10_pink.CL_from_OL_direct, opts1);
bodeplot(results_P15_pink.CL_from_OL_direct, opts1);
bodeplot(results_wc10_pink.CL_from_OL_direct, opts1);
bodeplot(results_wc15_pink.CL_from_OL_direct, opts1);
bodeplot(results_wc20_pink.CL_from_OL_direct, opts1);
legend('P7','P10','P15','wc10','wc15','wc20'); grid on;

figure; hold on; opts1.Title.String='CL';
% bodeplot(results_P7_sineSweep.CL, opts1);
% bodeplot(results_P15_sineSweep.CL, opts1);
% bodeplot(results_P5_pink.CL, opts1);
bodeplot(results_P7_pink.CL, opts1);
bodeplot(results_P10_pink.CL, opts1);
bodeplot(results_P15_pink.CL, opts1);
bodeplot(results_wc10_pink.CL, opts1);
bodeplot(results_wc15_pink.CL, opts1);
bodeplot(results_wc20_pink.CL, opts1);
legend('P7','P10','P15','wc10','wc15','wc20'); grid on;

figure; hold on; opts1.Title.String='CL from OL_{est nonLin}';
% bodeplot(results_P7_sineSweep.CL_from_OL_est_nonLin, opts1);
% bodeplot(results_P15_sineSweep.CL_from_OL_est_nonLin, opts1);
% bodeplot(results_P5_pink.CL_from_OL_est_nonLin, opts1);
bodeplot(results_P7_pink.CL_from_OL_est_nonLin, opts1);
bodeplot(results_P10_pink.CL_from_OL_est_nonLin, opts1);
bodeplot(results_P15_pink.CL_from_OL_est_nonLin, opts1);
bodeplot(results_wc10_pink.CL_from_OL_est_nonLin, opts1);
bodeplot(results_wc15_pink.CL_from_OL_est_nonLin, opts1);
bodeplot(results_wc20_pink.CL_from_OL_est_nonLin, opts1);
legend('P7','P10','P15','wc10','wc15','wc20'); grid on;





%%
% %fig9 = figure(9);ax9 = axes(fig9); hold(ax9, 'on');
% %bodeplot(G_closed,opts1);
% figure;hold on;
% bodeplot(G_PIDF_10Hz,opts1);    
% % bodeplot(G_PIDF_15Hz,opts1); 
% % bodeplot(G_PIDF_20Hz,opts1);    
% % bodeplot(Optimal_CL,opts1);
% % legend;
% bodeplot(results_wc10_pink.CL_from_OL_est_nonLin,opts1);    
% bodeplot(results_wc10_pink.CL_from_OL_direct,opts1);    
% bodeplot(results_wc10_pink.CL,opts1);   

% 
% %
% figure;hold on;
% bodeplot(G_PIDF_20Hz,opts1);    
% bodeplot(results_wc20_pink.CL_from_OL_est_nonLin,opts1);    
% bodeplot(results_wc20_pink.CL_from_OL_direct,opts1);    
% bodeplot(results_wc20_pink.CL,opts1);  


% %% Finding Response Spectre of Ground
% f_i=0.1; %freq inicial
% f_n=40;  %freq final
% n_points = 1e2;
% f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
% 
% [picos_ddx_ground , picos_x_ground] = ResponseSpectrum( t_vector , x_ref , ddx_ref, f_vector , 1);
% x_cl = lsim(d2d(G_closed,Ts) ,  x_ref , t_vector,'zoh');
% ddx_cl = secondDerivativeTime(x_cl , Ts);
% [picos_ddx_cl , picos_x_cl] = ResponseSpectrum( t_vector , x_cl , ddx_cl, f_vector , 1);
% 
% x_PIDF_5Hz = lsim(d2d(G_PIDF_20Hz,Ts) ,  x_ref , t_vector,'zoh');
% ddx_PIDF_5Hz = secondDerivativeTime(x_PIDF_5Hz , Ts);
% [picos_ddx_PIDF_5Hz , picos_x_PIDF_5Hz] = ResponseSpectrum( t_vector , x_PIDF_5Hz , ddx_PIDF_5Hz, f_vector , 1);
% 
% x_PIDF_10Hz = lsim(d2d(G_PIDF_10Hz,Ts) ,  x_ref , t_vector,'zoh');
% ddx_PIDF_10Hz = secondDerivativeTime(x_PIDF_10Hz , Ts);
% [picos_ddx_PIDF_10Hz , picos_x_PIDF_10Hz] = ResponseSpectrum( t_vector , x_PIDF_10Hz , ddx_PIDF_10Hz, f_vector , 1);
% 
% x_PIDF_15Hz = lsim(d2d(G_PIDF_15Hz,Ts) ,  x_ref , t_vector,'zoh');
% ddx_PIDF_15Hz = secondDerivativeTime(x_PIDF_15Hz , Ts);
% [picos_ddx_PIDF_15Hz , picos_x_PIDF_15Hz] = ResponseSpectrum( t_vector , x_PIDF_15Hz , ddx_PIDF_15Hz, f_vector , 1);
% 
% x_Optimal = lsim(d2d(Optimal_CL,Ts) ,  x_ref , t_vector,'zoh');
% ddx_Optimal = secondDerivativeTime(x_Optimal , Ts);
% [picos_ddx_Optimal , picos_x_Optimal] = ResponseSpectrum( t_vector , x_Optimal , ddx_Optimal, f_vector , 1);
% 
% baseFolder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\figures\';   % Base folder where you want to create the timestamped subfolder
% ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
% timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
% if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
%     mkdir(timeDir)
% end
% 
% fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);%ylim([0 ceil(max(picos_ddx_T_tuned(1:385,1))) ])
% subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
% color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4
% 
% figure(fig8); subplot(121); grid on; legend(); hold on;
% plot(f_vector, picos_ddx_ground,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
% plot(f_vector, picos_ddx_cl,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Standard -  MSE= %.2e',      mean((picos_ddx_ground-picos_ddx_cl).^2 )));% - Normal
% plot(f_vector, picos_ddx_PIDF_10Hz,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'PIDF 10Hz - MSE= %.2e',   mean((picos_ddx_ground-picos_ddx_PIDF_10Hz).^2 )));
% plot(f_vector, picos_ddx_PIDF_15Hz ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'PIDF 15Hz - MSE= %.2e', mean((picos_ddx_ground-picos_ddx_PIDF_15Hz).^2 )));
% plot(f_vector, picos_ddx_Optimal ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Optimal Control -  MSE= %.2e', mean((picos_ddx_ground-picos_ddx_Optimal).^2 )));
% 
% subplot(122); grid on;legend();hold on;
% plot(f_vector, picos_x_ground,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
% plot(f_vector, picos_x_cl,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Standard - MSE= %.2e',     mean((picos_x_ground-picos_x_cl).^2 )));%- Normal
% plot(f_vector, picos_x_PIDF_10Hz,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'PIDF 10Hz - MSE= %.2e',  mean((picos_x_ground-picos_x_PIDF_10Hz).^2 )));
% plot(f_vector, picos_x_PIDF_15Hz, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'PIDF 15Hz - MSE= %.2e', mean((picos_x_ground-picos_x_PIDF_15Hz).^2 )));
% plot(f_vector, picos_x_Optimal, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Optimal Control - MSE= %.2e', mean((picos_x_ground-picos_x_Optimal).^2 )));
% 
% set(fig8, 'WindowState', 'maximized');
% exportgraphics(fig8,fullfile(timeDir,'Response_Spectra_N.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');
% 
% %%
% figure; hold on; grid on;
% plot(t_vector, x_ref);
% plot(t_vector, x_cl);
% plot(t_vector, x_PIDF_15Hz);
% plot(t_vector , x_Optimal)