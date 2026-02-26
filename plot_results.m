%clear;clc;close all; 
setappdata(0, 'AutoStagger_LRDown_Last', []); set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\'; addpath(func_folder);
Ts = 0.005;Ts_fpga = 1/5000;
fir_np=100; np_CL=4; np_OL=4;

% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;%opts1.Ylim={[-40 10]};

% Control channel AI2 Displacement - 16 bit signed integer to mm conversion
a = 0.000485; b = -0.2;
bits2mm = @(bits) a*bits+b;
mm2bits = @(mm) (mm-b)/a;
clear a b

% % w_c = 10Hz
% Kp = 6.75;
% Ki = 68;
% Kd = .0712891;
% Tf = 0.00161743;
% pidf_wc10 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
% 
% % Data wc_15Hz
% Kp = 8.75
% Ki = 139
% Kd = 0.125977
% Tf = 0.00109863
% pidf_wc15 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
% 
% % w_c = 20Hz
% Kp = 8.562500;
% Ki = 75;
% Kd = 0.18457;
% Tf = 0.000167847;
% pidf_wc20 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)

%% input file - pink noise 40hz
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);

input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
clear x_drv_L_0  x_drv_V_0
%%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=15;
results_P15_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

OL_5000 = d2d(ss(results_P15_pink.OL_est_nonLin),Ts_fpga,'tustin');

%%
tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');

cutoff_frequency = 10; % Hz
PIDF  = pidtune(OL_5000,'PIDF',cutoff_frequency*2*pi,tuner_opts)
CL_PIDF_10Hz = d2d(feedback(PIDF*OL_5000, 1),Ts,'tustin');
cutoff_frequency = 15; % Hz
PIDF   = pidtune(OL_5000,'PIDF',cutoff_frequency*2*pi,tuner_opts)
CL_PIDF_15Hz = d2d(feedback(PIDF*OL_5000, 1),Ts,'tustin');
cutoff_frequency = 20; % Hz
PIDF  = pidtune(OL_5000,'PIDF',cutoff_frequency*2*pi,tuner_opts)
CL_PIDF_20Hz = d2d(feedback(PIDF*OL_5000, 1),Ts,'tustin');

%%
n_states=size(OL_5000.A,1);
OL_5000.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);

plant_aug = ss(OL_5000.A, OL_5000.B,[eye(n_states);OL_5000.C],0 , Ts_fpga);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [OL_5000.StateName ; {'y_xT'}];  % plant output

sumblk1 = sumblk('e = x_ref - y_xT'); % Compute the error signal: e = r - y

integrator = tf(1,[1 0], Ts_fpga); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error
%%
close all;
Q = diag([1e-9*ones(1,n_states),3e2]);%1e3*diag([zeros(1,n_states),1]);%blkdiag(eye(nx), eye(ny));
R = eye(size(OL_5000.B,2));
K_lqi = lqi(OL_5000, Q, R)% Design the LQI controller for the original system

K  = K_lqi(1:n_states);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = [ OL_5000.StateName ; {'xi'}];
controller.OutputName = {'i_sv'};
Optimal_CL = connect(plant_aug,  controller , integrator, sumblk1, 'x_ref','y_xT');

%%  Tolmezzo tgt
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\';
file = 'TolmezzoReducedScale.tgt'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', folder_1201); % load input drv
scale = 0.16;
x_tgt_T = scale*x_tgt_T; % convert to mm
ddx_tgt_T = scale*ddx_tgt_T; % convert to mm

% x_sim_T_10hz = lsim(CL_PIDF_10Hz ,  x_tgt_T , time_vector,'zoh');
% ddx_sim_T_10hz = secondDerivativeTime(x_sim_T_10hz , Ts);

x_sim_T_15hz = lsim(CL_PIDF_15Hz ,  x_tgt_T , time_vector,'zoh');
ddx_sim_T_15hz = secondDerivativeTime(x_sim_T_15hz , Ts);

% x_sim_T_20hz = lsim(CL_PIDF_20Hz ,  x_tgt_T , time_vector,'zoh');
% ddx_sim_T_20hz = secondDerivativeTime(x_sim_T_20hz , Ts);

x_sim_T_optimal = lsim(d2d(Optimal_CL,Ts) ,  x_tgt_T, time_vector,'zoh');
ddx_sim_T_optimal = secondDerivativeTime(x_sim_T_optimal , Ts);

% % Tolmezzo tune 10hz, 15hz , 20hz results
% folder_1201_tolmezzo ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\tolmezzo_scl0.16\';
% file = 'tune_cutoff_10Hz.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
% x_acq_T_10hz = x_acq_T;
% %sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% ddx_acq_T_10hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

% file = 'tune_cutoff_15hz.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
% x_acq_T_15hz= x_acq_T;
% %sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% ddx_acq_T_15hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

% file = 'tune_cutoff_20hz.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
% x_acq_T_20hz = x_acq_T;
% %sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% ddx_acq_T_20hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

%% Computing Response spectra of Adapted
folder_1201_tolmezzo_drv_acq ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\tolmezzo_drv&acq\';
% file = 'TolmezzoReducedScale_1.acq';
% file = 'TolmezzoReducedScale_2.acq';
file = 'TolmezzoReducedScale_4.acq'; 
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo_drv_acq );
ddx_acq_T_comp = secondDerivativeTime(x_acq_T , Ts);

%% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);

% Finding Response Spectre  of tuned pidf 10h 15hz 20hz in simulation
% [picos_ddx_sim_T_10hz, picos_x_sim_T_10hz] = ResponseSpectrum( time_vector , x_sim_T_10hz , ddx_sim_T_10hz, f_vector , 1);
[picos_ddx_sim_T_15hz, picos_x_sim_T_15hz] = ResponseSpectrum( time_vector , x_sim_T_15hz , ddx_sim_T_15hz, f_vector , 1);
% [picos_ddx_sim_T_20hz, picos_x_sim_T_20hz] = ResponseSpectrum( time_vector , x_sim_T_20hz , ddx_sim_T_20hz, f_vector , 1);

[picos_ddx_sim_T_optimal , picos_x_sim_T_optimal] = ResponseSpectrum( time_vector , x_sim_T_optimal , ddx_sim_T_optimal, f_vector , 1);

% Finding Response Spectre  of tuned pidf 10h 15hz 20hz with acquired data
% [picos_ddx_acq_T_10hz, picos_x_acq_T_10hz] = ResponseSpectrum( time_vector , x_acq_T_10hz , ddx_acq_T_10hz_comp, f_vector , 1);
% [picos_ddx_acq_T_15hz, picos_x_acq_T_15hz] = ResponseSpectrum( time_vector , x_acq_T_15hz , ddx_acq_T_15hz_comp, f_vector , 1);
% [picos_ddx_acq_T_20hz, picos_x_acq_T_20hz] = ResponseSpectrum( time_vector , x_acq_T_20hz , ddx_acq_T_20hz_comp, f_vector , 1);

%% Adapted driver response spectra
% [picos_ddx_acq_T_1  , picos_x_acq_T_1 ] = ResponseSpectrum( time_acq , x_acq_T, ddx_acq_T_comp, f_vector , 1);
% [picos_ddx_acq_T_acq_2  , picos_x_acq_T_acq_2 ] = ResponseSpectrum( time_vector , x_acq_T, ddx_acq_T_comp, f_vector , 1);
 [picos_ddx_acq_T_acq_4  , picos_x_acq_T_acq_4 ] = ResponseSpectrum( time_vector , x_acq_T, ddx_acq_T_comp, f_vector , 1);

%%
close;
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 10]);%ylim([0 ceil(max(picos_ddx_acq_T_20hz(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';color5='#77AC30';color6='#00fff7';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; hold on;
plot(f_vector, picos_ddx_tgt_T     ,'-' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
% plot(f_vector, picos_ddx_sim_T_10hz,'.-', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_10hz).^2 )));% - Normal
plot(f_vector, picos_ddx_sim_T_15hz,'.-', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_15hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_sim_T_20hz,'.-', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_20hz).^2 )));% - Normal
plot(f_vector, picos_ddx_sim_T_optimal,'.-', 'LineWidth' , 2, 'Color', color6,'DisplayName',sprintf( 'LQI (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_optimal).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_10hz,'--', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (acq)-  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_acq_T_10hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_15hz,'--', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (acq)-  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_acq_T_15hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_20hz,'--', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (acq)-  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_acq_T_20hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_1,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Driver 1 (acq) - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_acq_T_1).^2 )));
% plot(f_vector, picos_ddx_acq_T_acq_2,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Driver 2 (acq) - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_acq_T_acq_2).^2 )));
plot(f_vector, picos_ddx_acq_T_acq_4,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Driver 4 (acq) - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_acq_T_acq_4).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T     ,'-' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
% plot(f_vector, picos_x_sim_T_10hz,'.-', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_10hz).^2 )));% - Normal
plot(f_vector, picos_x_sim_T_15hz,'.-', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_15hz).^2 )));% - Normal
% plot(f_vector, picos_x_sim_T_20hz,'.-', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_20hz).^2 )));% - Normal
plot(f_vector, picos_x_sim_T_optimal,'.-', 'LineWidth' , 2, 'Color', color6,'DisplayName',sprintf( 'LQI (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_optimal).^2 )));% - Normal
% plot(f_vector, picos_x_acq_T_10hz,'--', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (acq)- MSE= %.2e',     mean((picos_x_tgt_T-picos_x_acq_T_10hz).^2 )));%- Normal
% plot(f_vector, picos_x_acq_T_15hz,'--', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (acq)- MSE= %.2e',     mean((picos_x_tgt_T-picos_x_acq_T_15hz).^2 )));%- Normal
% plot(f_vector, picos_x_acq_T_20hz,'--', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (acq)- MSE= %.2e',     mean((picos_x_tgt_T-picos_x_acq_T_20hz).^2 )));%- Normal
% plot(f_vector, picos_x_acq_T_1,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Driver 1 (acq) - MSE= %.2e', mean((picos_x_tgt_T-picos_x_acq_T_1).^2 )));
% plot(f_vector, picos_x_acq_T_acq_2,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Driver 2 (acq) - MSE= %.2e', mean((picos_x_tgt_T-picos_x_acq_T_acq_2).^2 )));
plot(f_vector, picos_x_acq_T_acq_4,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Driver 4 (acq) - MSE= %.2e', mean((picos_x_tgt_T-picos_x_acq_T_acq_4).^2 )));

fontsize(scale=1.8) ;

%%  laquila tgt
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\';
file = 'laquilaReducedScale.tgt'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', folder_1201); % load input drv
scale = 0.18;
x_tgt_T = scale*x_tgt_T; % convert to mm
ddx_tgt_T = scale*ddx_tgt_T; % convert to mm

% x_sim_T_10hz = lsim(CL_PIDF_10Hz ,  x_tgt_T , time_vector,'zoh');
% ddx_sim_T_10hz = secondDerivativeTime(x_sim_T_10hz , Ts);

x_sim_T_15hz = lsim(CL_PIDF_15Hz ,  x_tgt_T , time_vector,'zoh');
ddx_sim_T_15hz = secondDerivativeTime(x_sim_T_15hz , Ts);

% x_sim_T_20hz = lsim(CL_PIDF_20Hz ,  x_tgt_T , time_vector,'zoh');
% ddx_sim_T_20hz = secondDerivativeTime(x_sim_T_20hz , Ts);

x_sim_T_optimal = lsim(d2d(Optimal_CL,Ts) ,  x_tgt_T, time_vector,'zoh');
ddx_sim_T_optimal = secondDerivativeTime(x_sim_T_optimal , Ts);

% % laquila tune 10hz, 15hz , 20hz results
% folder_1201_tolmezzo ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\laquila_scl0.18\';
% file = 'tune_cutoff_10Hz.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
% x_acq_T_10hz = x_acq_T;
% %sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% ddx_acq_T_10hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

% file = 'tune_cutoff_15hz.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
% x_acq_T_15hz= x_acq_T;
% %sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% ddx_acq_T_15hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

% file = 'tune_cutoff_20hz.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
% x_acq_T_20hz = x_acq_T;
% %sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% ddx_acq_T_20hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

%% Computing Response spectra of Adapted
folder_1201_laquila_drv_acq ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\laquila_drv&acq\';

file = 'laquilaReducedScale_1.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_laquila_drv_acq );
x_acq_T_1=x_acq_T;
ddx_acq_T_1_comp = secondDerivativeTime(x_acq_T_1 , Ts);

% file = 'laquilaReducedScale_2.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_laquila_drv_acq );
% x_acq_T_2=x_acq_T;
% ddx_acq_T_2_comp = secondDerivativeTime(x_acq_T_2 , Ts);

% file = 'laquilaReducedScale_3.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_laquila_drv_acq );
% x_acq_T_3=x_acq_T;
% ddx_acq_T_3_comp = secondDerivativeTime(x_acq_T_3 , Ts);

% file = 'laquilaReducedScale_4.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_1201_laquila_drv_acq );
% x_acq_T_4=x_acq_T;
% ddx_acq_T_4_comp = secondDerivativeTime(x_acq_T_4 , Ts);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);

% Finding Response Spectre  of tuned pidf 10h 15hz 20hz in simulation
% [picos_ddx_sim_T_10hz, picos_x_sim_T_10hz] = ResponseSpectrum( time_vector , x_sim_T_10hz , ddx_sim_T_10hz, f_vector , 1);
[picos_ddx_sim_T_15hz, picos_x_sim_T_15hz] = ResponseSpectrum( time_vector , x_sim_T_15hz , ddx_sim_T_15hz, f_vector , 1);
% [picos_ddx_sim_T_20hz, picos_x_sim_T_20hz] = ResponseSpectrum( time_vector , x_sim_T_20hz , ddx_sim_T_20hz, f_vector , 1);

[picos_ddx_sim_T_optimal , picos_x_sim_T_optimal] = ResponseSpectrum( time_vector , x_sim_T_optimal , ddx_sim_T_optimal, f_vector , 1);

% Finding Response Spectre  of tuned pidf 10h 15hz 20hz
% [picos_ddx_acq_T_10hz, picos_x_acq_T_10hz] = ResponseSpectrum( time_vector , x_acq_T_10hz , ddx_acq_T_10hz_comp, f_vector , 1);
% [picos_ddx_acq_T_15hz, picos_x_acq_T_15hz] = ResponseSpectrum( time_vector , x_acq_T_15hz , ddx_acq_T_15hz_comp, f_vector , 1);
% [picos_ddx_acq_T_20hz, picos_x_acq_T_20hz] = ResponseSpectrum( time_vector , x_acq_T_20hz , ddx_acq_T_20hz_comp, f_vector , 1);

% Adapted driver response spectra
[picos_ddx_acq_T_1  , picos_x_acq_T_1 ] = ResponseSpectrum( time_acq , x_acq_T_1, ddx_acq_T_1_comp, f_vector , 1);
% [picos_ddx_acq_T_2  , picos_x_acq_T_2 ] = ResponseSpectrum( time_vector , x_acq_T_2, ddx_acq_T_2_comp, f_vector , 1);
% [picos_ddx_acq_T_3  , picos_x_acq_T_3 ] = ResponseSpectrum( time_vector , x_acq_T_3, ddx_acq_T_3_comp, f_vector , 1);
% [picos_ddx_acq_T_4  , picos_x_acq_T_4 ] = ResponseSpectrum( time_vector , x_acq_T_4, ddx_acq_T_4_comp, f_vector , 1);

%%
fig9 = figure(9);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 10]);%ylim([0 ceil(max(picos_ddx_acq_T_20hz(1:395,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);

figure(fig9); subplot(121); grid on; hold on;
plot(f_vector, picos_ddx_tgt_T     ,'-' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
% plot(f_vector, picos_ddx_sim_T_10hz,'.-', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_10hz).^2 )));% - Normal
plot(f_vector, picos_ddx_sim_T_15hz,'.-', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_15hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_sim_T_20hz,'.-', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_20hz).^2 )));% - Normal
plot(f_vector, picos_ddx_sim_T_optimal,'.-', 'LineWidth' , 2, 'Color', color6,'DisplayName',sprintf( 'LQI (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_T_optimal).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_10hz,'--', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (acq)-  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_acq_T_10hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_15hz,'--', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (acq)-  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_acq_T_15hz).^2 )));% - Normal
% plot(f_vector, picos_ddx_acq_T_20hz,'--', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (acq)-  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_acq_T_20hz).^2 )));% - Normal
plot(f_vector, picos_ddx_acq_T_1 ,'-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Driver 1 - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_acq_T_1).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T     ,'-' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
% plot(f_vector, picos_x_sim_T_10hz,'.-', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_10hz).^2 )));% - Normal
plot(f_vector, picos_x_sim_T_15hz,'.-', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_15hz).^2 )));% - Normal
% plot(f_vector, picos_x_sim_T_20hz,'.-', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_20hz).^2 )));% - Normal
plot(f_vector, picos_x_sim_T_optimal,'.-', 'LineWidth' , 2, 'Color', color6,'DisplayName',sprintf( 'LQI (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_T_optimal).^2 )));% - Normal
% plot(f_vector, picos_x_acq_T_10hz,'--', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=10Hz (acq)- MSE= %.2e',     mean((picos_x_tgt_T-picos_x_acq_T_10hz).^2 )));%- Normal
% plot(f_vector, picos_x_acq_T_15hz,'--', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'PIDF ω_c=15Hz (acq)- MSE= %.2e',     mean((picos_x_tgt_T-picos_x_acq_T_15hz).^2 )));%- Normal
% plot(f_vector, picos_x_acq_T_20hz,'--', 'LineWidth' , 2, 'Color', color5,'DisplayName',sprintf( 'PIDF ω_c=20Hz (acq)- MSE= %.2e',     mean((picos_x_tgt_T-picos_x_acq_T_20hz).^2 )));%- Normal
plot(f_vector, picos_x_acq_T_1, '-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Driver 1 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_acq_T_1).^2 )));

fontsize(scale=1.8) ;

%%
baseFolder = folder_1201;   % Base folder where you want to create the timestamped subfolder
ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
    mkdir(timeDir)
end

set(fig8, 'WindowState', 'maximized');
exportgraphics(fig8,fullfile(timeDir,'Response_Spectra_Tolmezzo.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');

set(fig9, 'WindowState', 'maximized');
exportgraphics(fig9,fullfile(timeDir,'Response_Spectra_Laquila.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');