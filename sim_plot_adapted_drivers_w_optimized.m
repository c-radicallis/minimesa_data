clear;clc;close all;
addpath ('C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model' , ...
    'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\',...
    'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\functions_matlab')
%
launch_Adapt =0; % Set to 1 to lauch Adapt.exe
return_on = 1; % Set to 1 for execution to stop before adapting drivers, or set to 0 if the adapted drivers have already been generated

opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;
opts1.PhaseVisible='off'; opts1.YLim={[-30 10]};
 Ts = 0.005;

%% Load Plant Model
load('optimized_benchmark_results\OL200_from_P10.mat')
% % input file - pink noise 40hz
% input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
% file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
% x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
% clear x_drv_L_0  x_drv_V_0
% %  Data P10
% folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
% file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
% x_acq_T = x_acq_T*1e3;
% a = 0.000485; b = -0.2; bits2mm = @(bits) a*bits+b; mm2bits = @(mm) (mm-b)/a; clear a b; % Control channel AI2 Displacement - 16 bit signed integer to mm conversion
% sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% Kp=10; fir_np=100; np_CL=4; np_OL=4;
% id_results = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);
% OL_200 = ss(id_results.OL_est_nonLin)

%% Create augmented state space system
n_states = size(OL_200.A,1); % Create augmented state space model
OL_200.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);
plant_aug = ss(OL_200.A, OL_200.B,[eye(n_states);OL_200.C],[zeros(n_states,1); OL_200.D] , Ts);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [OL_200.StateName ; {'y_xT'}];  % plant output
sumblk1 = sumblk('e = x_ref - y_xT'); % Compute the error signal: e = r - y
integrator = tf(1,[1 -1], Ts);  integrator.InputName = {'e'};  integrator.OutputName = {'xi'};  % The integrator integrates the tracking error. % error: e = r - y % integrated error

%% Load target
folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\PRJ_Elcentro\';
target = 'Elcentro.tgt'; 
LTF_to_TXT_then_load(target,'InputFolder', folder)

disp_limit=4e-3; max_tgt=max([x_tgt_T ; x_tgt_L]);
if max_tgt>disp_limit
    scale = round(disp_limit/max_tgt , 2)
    x_tgt_T   = scale*x_tgt_T; x_tgt_L   = scale*x_tgt_L; ddx_tgt_T = scale*ddx_tgt_T; ddx_tgt_L = scale*ddx_tgt_L;
end
max_abs_x_tgt_T = max( abs(x_tgt_T ))
max_abs_x_tgt_L = max( abs(x_tgt_L ))

%% ── Run the optimisation of Q
Q1_0 = 4.5183e+14; Q2_0 = 3.4447e-01; Q3_0 = 4.0431e+03; Q4_0 = 5.1600e-02; Qi_0 =3.1503e+15;  % for Kp = 10 % Initial guess (log-space) — edit to reflect your engineering intuition
% Q1_0 = 1;  Q2_0 = 1;  Q3_0 = 1;  Q4_0 = 1;  Qi_0 = 10;

log_q0 = log([Q1_0, Q2_0, Q3_0, Q4_0, Qi_0]);
outputFcn = @(~, ov, state) recordAndStop(ov, state);
opts_opt = optimset('Display',     'iter', ...
                    'TolX',        1e-5,  ...
                    'TolFun',      1e-5,  ...
                    'MaxFunEvals', 1e12,   ...
                    'MaxIter',     1e12,   ...
                    'OutputFcn',   outputFcn);

picos_ddx_tgt_T_ForCost = ResponseSpectrumForCost(  ddx_tgt_T );
objFun = @(log_q) AccelSpectraCost(log_q, OL_200, plant_aug, integrator, sumblk1,   n_states, picos_ddx_tgt_T_ForCost , ddx_tgt_T , time_vector );

fprintf('=== Starting Q optimisation ===\n');
[log_q_best, J_best] = fminsearch(objFun, log_q0, opts_opt);

% ── Recover & display best weights ───────────────────────────────────────
q_best = exp(log_q_best); Q1_best = q_best(1); Q2_best = q_best(2); Q3_best = q_best(3); Q4_best = q_best(4); Qi_best = q_best(5);
fprintf('\n=== Optimisation complete ===\n');
fprintf('Best Q weights:\n');
fprintf('  Q1 = %.4e\n  Q2 = %.4e\n  Q3 = %.4e\n  Q4 = %.4e\n  Qi = %.4e\n', Q1_best, Q2_best, Q3_best, Q4_best, Qi_best);

Q_best   = diag(q_best);
R        = 1;
K_lqi_best = lqi(OL_200, Q_best, R)
K        = K_lqi_best(1:n_states);
Ki       = K_lqi_best(end);

controller_best = ss([], [], [], -[K, Ki]);
controller_best.InputName  = [OL_200.StateName; {'xi'}];
controller_best.OutputName = {'i_sv'};

CL_LQI = connect(plant_aug, controller_best, integrator, sumblk1, 'x_ref', 'y_xT');
figure; hold on; bodeplot(CL_PIDF, CL_LQI , opts1); grid on;legend;

x_T_LQI = lsim(CL_LQI ,  x_tgt_T , time_vector,'zoh');
ddx_T_LQI = secondDerivativeTime(x_T_LQI , Ts);
x_L_LQI = lsim(CL_LQI ,  x_tgt_L , time_vector,'zoh');
ddx_L_LQI = secondDerivativeTime(x_L_LQI , Ts);

%% Tuning PIDF  and running simulations
tuner_opts = pidtuneOptions('DesignFocus','reference-tracking'); % Tune PIDF
cutoff_frequency = 7; % Hz
PIDF   = pidtune(OL_200,'PIDF',cutoff_frequency*2*pi,tuner_opts)
CL_PIDF = feedback(PIDF*OL_200, 1);
x_T_tuned = lsim(CL_PIDF ,  x_tgt_T , time_vector,'zoh');
ddx_T_tuned = secondDerivativeTime(x_T_tuned , Ts);
x_L_tuned = lsim(CL_PIDF ,  x_tgt_L , time_vector,'zoh');
ddx_L_tuned = secondDerivativeTime(x_L_tuned , Ts);

%% Lauch Adapt.exe % note the empty quotes "" are the window title placeholder
if launch_Adapt
    cmd = sprintf('start "" "%s"', fullfile('C:','Users','afons','OneDrive - Universidade de Lisboa','Controlo de Plataforma Sismica','uniaxial_table_model','Adapting_Driver_Signal','Adapt.exe.lnk'));
    system(cmd);
    fprintf("Launched Adapt.exe, continuing script...\n \n ");
end
fprintf("\n \n Go to Adapt.exe and generate driver 0 (Click 'Adapt Init' button) \n \n ")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 0
name = target(1 : end-4); %#ok<UNRCH>
LTF_to_TXT_then_load( [ name, '_0.DRV' ] ,'InputFolder',folder)
x_T_acq_0 = lsim(G_xT_xref ,  x_drv_T_0 , time_vector,'zoh');
ddx_T_acq_0 = secondDerivativeTime(x_T_acq_0 , Ts);
% writeTXT_then_LTF(time_vector,x_T_acq_0,ddx_T_acq_0,folder,[ name, '_0.ACQ.txt' ]);
x_L_acq_0 = lsim(G_xT_xref ,  x_drv_L_0 , time_vector,'zoh');
ddx_L_acq_0 = secondDerivativeTime(x_L_acq_0 , Ts);
writeTXT_then_LTF(time_vector,[x_T_acq_0,x_L_acq_0],[ddx_T_acq_0,ddx_L_acq_0],folder,[ name, '_0.ACQ.txt' ]);
fprintf("\n \n Go to Adapt.exe and generate driver 1 (Click 'Process' button)\n \n")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 1
LTF_to_TXT_then_load( [ name, '_1.DRV' ] ,'InputFolder',folder)
x_T_acq_1 = lsim(G_xT_xref ,  x_drv_T_1 , time_vector,'zoh');
ddx_T_acq_1 = secondDerivativeTime(x_T_acq_1 , Ts);
%writeTXT_then_LTF(time_vector,x_T_acq_1,ddx_T_acq_1,folder,[ name, '_1.ACQ.txt' ]);
x_L_acq_1 = lsim(G_xT_xref ,  x_drv_L_1 , time_vector,'zoh');
ddx_L_acq_1 = secondDerivativeTime(x_L_acq_1 , Ts);
writeTXT_then_LTF(time_vector,[x_T_acq_1,x_L_acq_1],[ddx_T_acq_1,ddx_L_acq_1],folder, [ name, '_1.ACQ.txt' ]); 
fprintf("\n \n Go to Adapt.exe and generate driver 2 (Click 'Next Iteration' and then 'Process' button) \n \n")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 2
LTF_to_TXT_then_load( [ name, '_2.DRV' ] ,'InputFolder',folder)
x_T_acq_2 = lsim(G_xT_xref ,  x_drv_T_2 , time_vector,'zoh');
ddx_T_acq_2 = secondDerivativeTime(x_T_acq_2 , Ts);
% writeTXT_then_LTF(time_vector,x_T_acq_2,ddx_T_acq_2,folder,[ name, '_2.ACQ.txt' ]);
x_L_acq_2 = lsim(G_xT_xref ,  x_drv_L_2 , time_vector,'zoh');
ddx_L_acq_2 = secondDerivativeTime(x_L_acq_2 , Ts);
writeTXT_then_LTF(time_vector,[x_T_acq_2,x_L_acq_2],[ddx_T_acq_2,ddx_L_acq_2],folder, [ name, '_2.ACQ.txt' ]); 

%% Response Spectra settings
f_i=0.1; %freq inicial
f_n=20;  %freq final
n_points = 5e2;
f_vector_accel = logspace( log10(f_i) , log10(f_n) , n_points);
f_vector_disp = f_vector_accel(1:200);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum(f_vector_accel, ddx_tgt_T, x_tgt_T, f_vector_disp);
[picos_ddx_tgt_L , picos_x_tgt_L] = ResponseSpectrum(f_vector_accel, ddx_tgt_L, x_tgt_L, f_vector_disp);

% Response Spectre  of Optimal
[picos_ddx_T_tuned , picos_x_T_tuned] = ResponseSpectrum( f_vector_accel, ddx_T_tuned, x_T_tuned , f_vector_disp);
[picos_ddx_L_tuned , picos_x_L_tuned] = ResponseSpectrum( f_vector_accel, ddx_L_tuned, x_L_tuned , f_vector_disp);

% Response Spectre  of Optimal
[picos_ddx_T_LQI , picos_x_T_LQI] = ResponseSpectrum( f_vector_accel, ddx_T_LQI, x_T_LQI , f_vector_disp);
[picos_ddx_L_LQI , picos_x_L_LQI] = ResponseSpectrum( f_vector_accel, ddx_L_LQI, x_L_LQI , f_vector_disp);

%% Computing Response spectra of Adapted
[picos_ddx_T_acq_0  , picos_x_T_acq_0 ] = ResponseSpectrum( f_vector_accel, ddx_T_acq_0 , x_T_acq_0, f_vector_disp);
[picos_ddx_T_acq_1  , picos_x_T_acq_1 ] = ResponseSpectrum( f_vector_accel, ddx_T_acq_1 , x_T_acq_1, f_vector_disp);
[picos_ddx_T_acq_2  , picos_x_T_acq_2 ] = ResponseSpectrum( f_vector_accel, ddx_T_acq_2 , x_T_acq_2, f_vector_disp);
[picos_ddx_L_acq_0  , picos_x_L_acq_0 ] = ResponseSpectrum( f_vector_accel, ddx_L_acq_0 , x_L_acq_0, f_vector_disp);
[picos_ddx_L_acq_1  , picos_x_L_acq_1 ] = ResponseSpectrum( f_vector_accel, ddx_L_acq_1 , x_L_acq_1, f_vector_disp);
[picos_ddx_L_acq_2  , picos_x_L_acq_2 ] = ResponseSpectrum( f_vector_accel, ddx_L_acq_2 , x_L_acq_2, f_vector_disp);

%%
fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);ylim([0 ceil(max(picos_ddx_T_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_T_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_T_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',   mean((picos_ddx_tgt_T-picos_ddx_T_LQI).^2 )));
% plot(f_vector, picos_ddx_T_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_0).^2 )));
% plot(f_vector, picos_ddx_T_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_1).^2 )));
% plot(f_vector, picos_ddx_T_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_T_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned).^2 )));%- Normal
plot(f_vector, picos_x_T_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',  mean((picos_x_tgt_T-picos_x_T_LQI).^2 )));
% plot(f_vector, picos_x_T_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_0).^2 )));
% plot(f_vector, picos_x_T_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_1).^2 )));
% plot(f_vector, picos_x_T_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_2).^2 )));

% Create Figures - Longitudinal
fig9 = figure(9);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Parallel');xlim([1 20]);ylim([0 ceil(max(picos_ddx_L_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Parallel');xlim([0.1 5]);

figure(fig9); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_L,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_L_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e',      mean((picos_ddx_tgt_L-picos_ddx_L_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_L_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',   mean((picos_ddx_tgt_L-picos_ddx_L_LQI).^2 )));
% plot(f_vector, picos_ddx_L_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_0).^2 )));
% plot(f_vector, picos_ddx_L_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_1).^2 )));
% plot(f_vector, picos_ddx_L_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_L,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_L_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e',     mean((picos_x_tgt_L-picos_x_L_tuned).^2 )));%- Normal
plot(f_vector, picos_x_L_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',  mean((picos_x_tgt_L-picos_x_L_LQI).^2 )));
% plot(f_vector, picos_x_L_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_0).^2 )));
% plot(f_vector, picos_x_L_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_1).^2 )));
% plot(f_vector, picos_x_L_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_2).^2 )));

%% Simulation using updated driver 3
% LTF_to_TXT_then_load( [ name, '_3.DRV' ] ,'InputFolder',folder)
% x_T_acq_3 = lsim(G_xT_xref ,  x_drv_T_3 , time_vector,'zoh');
% ddx_T_acq_3 = secondDerivativeTime(x_T_acq_3 , Ts);
% % writeTXT_then_LTF(time_vector,x_T_acq_3,ddx_T_acq_3,folder,[ name, '_3.ACQ.txt' ]);
% x_L_acq_3 = lsim(G_xT_xref ,  x_drv_L_3 , time_vector,'zoh');
% ddx_L_acq_3 = secondDerivativeTime(x_L_acq_3 , Ts);
% writeTXT_then_LTF(time_vector,[x_T_acq_3,x_L_acq_3],[ddx_T_acq_3,ddx_L_acq_3],folder, [ name, '_3.ACQ.txt' ]); 
% 
% [picos_ddx_T_acq_3  , picos_x_T_acq_3 ] = ResponseSpectrum( time_vector , x_T_acq_3 , ddx_T_acq_3, f_vector , 1);
% % Create Figures - Trnaversal
% figure(fig8); subplot(121);
% plot(f_vector, picos_ddx_T_acq_3 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 3 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_3).^2 )));
% subplot(122);
% plot(f_vector, picos_x_T_acq_3, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 3 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_3).^2 )));
% 
% [picos_ddx_L_acq_3  , picos_x_L_acq_3 ] = ResponseSpectrum( time_vector , x_L_acq_3 , ddx_L_acq_3, f_vector , 1);
% % Create Figures - Longitudinal
% figure(fig9); subplot(121);
% plot(f_vector, picos_ddx_L_acq_3 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 3 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_3).^2 )));
% subplot(122);
% plot(f_vector, picos_x_L_acq_3, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 3 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_3).^2 )));


%%

baseFolder = folder;   % Base folder where you want to create the timestamped subfolder
ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
    mkdir(timeDir)
end

set(fig8, 'WindowState', 'maximized');
exportgraphics(fig8,fullfile(timeDir,'Response_Spectra_N.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');

set(fig9, 'WindowState', 'maximized');
exportgraphics(fig9,fullfile(timeDir,'Response_Spectra_P.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');