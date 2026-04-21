%% ── Q-weight optimisation via fminsearch ─────────────────────────────────
% Weights are searched in log-space so the optimiser can roam over many
% orders of magnitude without going negative.

clear;clc;close all; 

%%
setappdata(0, 'AutoStagger_LRDown_Last', []); set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\functions_matlab'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\'; addpath(func_folder);

Ts = 0.005;fir_np=100; np_CL=4; np_OL=4;

opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;
opts1.PhaseVisible='off'; opts1.YLim={[-30 10]};

a = 0.000485; b = -0.2; bits2mm = @(bits) a*bits+b; mm2bits = @(mm) (mm-b)/a; clear a b; % Control channel AI2 Displacement - 16 bit signed integer to mm conversion

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
clear x_drv_L_0  x_drv_V_0
%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=15;
results_P15_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);
OL_200 = ss(results_P15_pink.OL_est_nonLin)

tuner_opts = pidtuneOptions('DesignFocus','reference-tracking'); % Tune PIDF
cutoff_frequency = 15; % Hz
PIDF   = pidtune(OL_200,'PIDF',cutoff_frequency*2*pi,tuner_opts)
CL_PIDF_15Hz = feedback(PIDF*OL_200, 1);

n_states = size(OL_200.A,1); % Create augmented state space model
OL_200.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);

plant_aug = ss(OL_200.A, OL_200.B,[eye(n_states);OL_200.C],[zeros(n_states,1); OL_200.D] , Ts);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [OL_200.StateName ; {'y_xT'}];  % plant output

sumblk1 = sumblk('e = x_ref - y_xT'); % Compute the error signal: e = r - y
integrator = tf(1,[1 -1], Ts);  integrator.InputName = {'e'};  integrator.OutputName = {'xi'};  % The integrator integrates the tracking error. % error: e = r - y % integrated error

%%  load Tolmezzo tgt
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\';
file = 'TolmezzoReducedScale.tgt'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', folder_1201); % load input drv
scale = 0.16;
x_tgt_T = scale*x_tgt_T; 
ddx_tgt_T = scale*ddx_tgt_T; 

x_sim_PIDF = lsim(CL_PIDF_15Hz, x_tgt_T, time_vector, 'zoh');     x_MSE_PIDF = mean((x_sim_PIDF - x_tgt_T).^2);
ddx_sim_PIDF = lsim(CL_PIDF_15Hz, ddx_tgt_T, time_vector, 'zoh'); ddx_MSE_PIDF = mean(( ddx_sim_PIDF -  ddx_tgt_T).^2);

%% ── Run the optimisation ─────────────────────────────────────────────────

% Initial guess (log-space) — edit to reflect your engineering intuition
 % Q1_0 = 1e-1;  Q2_0 = 1e-1;  Q3_0 = 1e-1;  Q4_0 = 1e-1;  Qi_0 = 1e1;
% Q1_0 = eps;  Q2_0 = eps;  Q3_0 = eps;  Q4_0 = eps;  Qi_0 = 200;
Q1_0 = 1;  Q2_0 = 1;  Q3_0 = 1;  Q4_0 = 1;  Qi_0 = 1;
log_q0 = log([Q1_0, Q2_0, Q3_0, Q4_0, Qi_0]);

outputFcn = @(~, ov, state) recordAndStop(ov, state);

opts_opt = optimset('Display',     'iter', ...
                    'TolX',        1e-15,  ...
                    'TolFun',      1e-15,  ...
                    'MaxFunEvals', 1e12,   ...
                    'MaxIter',     1e12,   ...
                    'OutputFcn',   outputFcn);

% Wrap objective so fminsearch only sees log_q
% objFun = @(log_q) trackingCost(log_q, OL_200, plant_aug, integrator, sumblk1, x_tgt_T, time_vector, n_states);
objFun = @(log_q) AccelTrackingCost(log_q, OL_200, plant_aug, integrator, sumblk1, ddx_tgt_T, time_vector, n_states);

fprintf('=== Starting Q optimisation ===\n');
[log_q_best, J_best] = fminsearch(objFun, log_q0, opts_opt);

% ── Recover & display best weights ───────────────────────────────────────
q_best = exp(log_q_best); Q1_best = q_best(1); Q2_best = q_best(2); Q3_best = q_best(3); Q4_best = q_best(4); Qi_best = q_best(5);
fprintf('\n=== Optimisation complete ===\n');
% fprintf('PIDF ddx_MSE : %.6e\n', ddx_MSE_PIDF);
% fprintf('Best LQI ddx_MSE : %.6f\n', J_best);
fprintf('Best Q weights:\n');
fprintf('  Q1 = %.4e\n  Q2 = %.4e\n  Q3 = %.4e\n  Q4 = %.4e\n  Qi = %.4e\n', Q1_best, Q2_best, Q3_best, Q4_best, Qi_best);

figure;
semilogy(opt_hist(:,1), opt_hist(:,2), 'b-o', 'MarkerSize', 3, 'LineWidth', 1); hold on;
yline(ddx_MSE_PIDF, 'r--', 'LineWidth', 1.5, 'Label', 'PIDF baseline');
yline(J_best,   'g--', 'LineWidth', 1.5, 'Label', 'Best LQI');
xlabel('Iteration');  ylabel('MSE (log scale)');
title('Optimisation convergence');  grid on;

% ── Rebuild and plot the best closed-loop system ─────────────────────────
Q_best   = diag(q_best);
R        = 1;
K_lqi_best    = lqi(OL_200, Q_best, R)
K        = K_lqi_best(1:n_states);
Ki       = K_lqi_best(end);

controller_best = ss([], [], [], -[K, Ki]);
controller_best.InputName  = [OL_200.StateName; {'xi'}];
controller_best.OutputName = {'i_sv'};

Optimal_CL_best = connect(plant_aug, controller_best, integrator, sumblk1, 'x_ref', 'y_xT');

x_sim_best = lsim(Optimal_CL_best, x_tgt_T, time_vector, 'zoh');    x_MSE_best = mean((x_sim_best - x_tgt_T).^2);
ddx_sim_best = lsim(Optimal_CL_best, ddx_tgt_T, time_vector, 'zoh');ddx_MSE_best = mean(( ddx_sim_best -  ddx_tgt_T).^2);

%% Manual tunning
Q_manual = diag([ones(1,n_states) , 1e1 ]) %eps*eye(1,n_states),
R = eye(size(OL_200.B,2));
K_lqi_manual = lqi(OL_200, Q_manual, R)% Design the LQI controller for the original system

K  = K_lqi_manual(1:n_states);      % state feedback gains
Ki = K_lqi_manual(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); controller.InputName = [ OL_200.StateName ; {'xi'}]; controller.OutputName = {'i_sv'}; %   u = -[K  Ki] * [x; xi]
Optimal_CL_manual = connect(plant_aug,  controller , integrator, sumblk1, 'x_ref','y_xT');
Optimal_CL_manual_stable = isstable(Optimal_CL_manual);
x_sim_manual = lsim(Optimal_CL_manual, x_tgt_T, time_vector, 'zoh'); x_MSE_manual=mean((x_sim_manual - x_tgt_T).^2);
ddx_sim_manual = lsim(Optimal_CL_manual, ddx_tgt_T, time_vector, 'zoh'); ddx_MSE_manual=mean((ddx_sim_manual - ddx_tgt_T).^2);

%%
close all;
figure;hold on;
plot(time_vector, x_tgt_T,  'g','LineWidth', 1);  set(gca, 'ColorOrderIndex', 1);  % reset counter
plot(time_vector, x_sim_PIDF, '-',  'LineWidth', 1);
plot(time_vector, x_sim_best,'-',  'LineWidth', 1);
plot(time_vector, x_sim_manual, '-',  'LineWidth', 1);
xlabel('Time (s)');  ylabel('x_T');
legend('Target', ...
       sprintf('PIDF (x_MSE = %.2e)',        x_MSE_PIDF), ...
       sprintf('Optimised LQI (x_MSE = %.2e)', x_MSE_best), ...
       sprintf('Manual LQI (x_MSE = %.2e)',  x_MSE_manual));
grid on;
xlim([2.5 5]); ylim('auto');

figure; hold on;
bodeplot(CL_PIDF_15Hz, opts1);
bodeplot(Optimal_CL_best, opts1);
bodeplot(Optimal_CL_manual, opts1);
title('Bode - Optimised closed-loop'); grid on;legend;

%% Response Spectra

% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);
[picos_ddx_sim_PIDF, picos_x_sim_PIDF] = ResponseSpectrum( time_vector , x_sim_PIDF , ddx_sim_PIDF, f_vector , 1);
[picos_ddx_sim_best , picos_x_sim_best] = ResponseSpectrum( time_vector , x_sim_best , ddx_sim_best, f_vector , 1);
[picos_ddx_sim_manual , picos_x_sim_manual] = ResponseSpectrum( time_vector , x_sim_manual , ddx_sim_manual, f_vector , 1);

%%
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 30]);%ylim([0 ceil(max(picos_ddx_acq_T_20hz(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 10]);
% color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';color5='#77AC30';color6='#00fff7';% Define colors for lines 1/3 and 2/4
RGB =  get(groot,"FactoryAxesColorOrder"); H = compose("#%02X%02X%02X",round(RGB*255)); color1 = 'g';color2 = H(1) ;color3 = H(2) ; color4 = H(3);

figure(fig8); subplot(121); grid on; hold on;legend()
plot(f_vector, picos_ddx_tgt_T     ,'-' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_sim_PIDF,'.-', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=15Hz (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_PIDF).^2 )));% - Normal
plot(f_vector, picos_ddx_sim_best,'.-', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'Optimised LQI (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_best).^2 )));% - Normal
plot(f_vector, picos_ddx_sim_manual,'.-', 'LineWidth' , 2, 'Color', color4,'DisplayName',sprintf( 'Manual LQI (sim) -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_sim_manual).^2 )));% - Normal

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T     ,'-' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_sim_PIDF,'.-', 'LineWidth' , 2, 'Color', color2,'DisplayName',sprintf( 'PIDF ω_c=15Hz (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_PIDF).^2 )));% - Normal
plot(f_vector, picos_x_sim_best,'.-', 'LineWidth' , 2, 'Color', color3,'DisplayName',sprintf( 'Optimised LQI (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_best).^2 )));% - Normal
plot(f_vector, picos_x_sim_manual,'.-', 'LineWidth' , 2, 'Color', color4,'DisplayName',sprintf( 'Manual LQI (sim) -  MSE= %.2e',      mean((picos_x_tgt_T-picos_x_sim_manual).^2 )));% - Normal

fontsize(scale=1.8);set(fig8, 'WindowState', 'maximized');

%%
% ddx_sim = lsim(Optimal_CL_best, ddx_tgt_T, time_vector, 'zoh');
% ddx_sim_foh = lsim(Optimal_CL_best, ddx_tgt_T, time_vector, 'foh');
% figure; hold on, grid on;
% plot(time_vector,ddx_sim);
% plot(time_vector,ddx_sim_foh);
% plot(time_vector,ddx_sim_best);
% legend('from lsim zoh', 'from lsim foh','from secondDerivativeTime')