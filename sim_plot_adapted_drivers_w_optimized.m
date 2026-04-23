clear;clc;close all;
addpath ('C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model' , ...
    'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\')
%
launch_Adapt =0; % Set to 1 to lauch Adapt.exe
return_on = 1; % Set to 1 for execution to stop before adapting drivers, or set to 0 if the adapted drivers have already been generated

%% Load Plant Model
% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
clear x_drv_L_0  x_drv_V_0
%  Data P10
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
a = 0.000485; b = -0.2; bits2mm = @(bits) a*bits+b; mm2bits = @(mm) (mm-b)/a; clear a b; % Control channel AI2 Displacement - 16 bit signed integer to mm conversion
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=10; Ts = 0.005; fir_np=100; np_CL=4; np_OL=4;
id_results = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);
OL_200 = ss(id_results.OL_est_nonLin)

%% Tuning PIDF
tuner_opts = pidtuneOptions('DesignFocus','reference-tracking'); % Tune PIDF
cutoff_frequency = 7; % Hz
PIDF   = pidtune(OL_200,'PIDF',cutoff_frequency*2*pi,tuner_opts)
CL_PIDF = feedback(PIDF*OL_200, 1);

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

limit=4e-3; max_tgt=max([x_tgt_T ; x_tgt_L]);
if max_tgt>limit
    scale = round(limit/max_tgt , 2)
    x_tgt_T   = scale*x_tgt_T; x_tgt_L   = scale*x_tgt_L; ddx_tgt_T = scale*ddx_tgt_T; ddx_tgt_L = scale*ddx_tgt_L;
end
max_abs_x_tgt_T = max( abs(x_tgt_T ))
max_abs_x_tgt_L = max( abs(x_tgt_L ))
% Compute target Response Spectra and define settings
f_i=0.1; %freq inicial
f_n=20;  %freq final
n_points = 1e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);



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
ddx_T_acq_0 = secondDerivativeTime(x_T_acq_0 , t_step);
% writeTXT_then_LTF(time_vector,x_T_acq_0,ddx_T_acq_0,folder,[ name, '_0.ACQ.txt' ]);
x_L_acq_0 = lsim(G_xT_xref ,  x_drv_L_0 , time_vector,'zoh');
ddx_L_acq_0 = secondDerivativeTime(x_L_acq_0 , t_step);
writeTXT_then_LTF(time_vector,[x_T_acq_0,x_L_acq_0],[ddx_T_acq_0,ddx_L_acq_0],folder,[ name, '_0.ACQ.txt' ]);
fprintf("\n \n Go to Adapt.exe and generate driver 1 (Click 'Process' button)\n \n")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 1
LTF_to_TXT_then_load( [ name, '_1.DRV' ] ,'InputFolder',folder)
x_T_acq_1 = lsim(G_xT_xref ,  x_drv_T_1 , time_vector,'zoh');
ddx_T_acq_1 = secondDerivativeTime(x_T_acq_1 , t_step);
%writeTXT_then_LTF(time_vector,x_T_acq_1,ddx_T_acq_1,folder,[ name, '_1.ACQ.txt' ]);
x_L_acq_1 = lsim(G_xT_xref ,  x_drv_L_1 , time_vector,'zoh');
ddx_L_acq_1 = secondDerivativeTime(x_L_acq_1 , t_step);
writeTXT_then_LTF(time_vector,[x_T_acq_1,x_L_acq_1],[ddx_T_acq_1,ddx_L_acq_1],folder, [ name, '_1.ACQ.txt' ]); 
fprintf("\n \n Go to Adapt.exe and generate driver 2 (Click 'Next Iteration' and then 'Process' button) \n \n")
if return_on
    return;
end   % execution stops here; lines below wonnt run

%% Simulation using updated driver 2
LTF_to_TXT_then_load( [ name, '_2.DRV' ] ,'InputFolder',folder)
x_T_acq_2 = lsim(G_xT_xref ,  x_drv_T_2 , time_vector,'zoh');
ddx_T_acq_2 = secondDerivativeTime(x_T_acq_2 , t_step);
% writeTXT_then_LTF(time_vector,x_T_acq_2,ddx_T_acq_2,folder,[ name, '_2.ACQ.txt' ]);
x_L_acq_2 = lsim(G_xT_xref ,  x_drv_L_2 , time_vector,'zoh');
ddx_L_acq_2 = secondDerivativeTime(x_L_acq_2 , t_step);
writeTXT_then_LTF(time_vector,[x_T_acq_2,x_L_acq_2],[ddx_T_acq_2,ddx_L_acq_2],folder, [ name, '_2.ACQ.txt' ]); 

%% Create Figures - Transversal
close all;

% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);
[picos_ddx_tgt_L , picos_x_tgt_L] = ResponseSpectrum( time_vector , x_tgt_L , ddx_tgt_L, f_vector , 1);

% Response Spectre  of Optimal
[picos_ddx_T_tuned , picos_x_T_tuned] = ResponseSpectrum( time_vector , x_T_tuned , ddx_T_tuned, f_vector , 1);
[picos_ddx_L_tuned , picos_x_L_tuned] = ResponseSpectrum( time_vector , x_L_tuned , ddx_L_tuned, f_vector , 1);

% Response Spectre  of Optimal
[picos_ddx_T_LQI , picos_x_T_LQI] = ResponseSpectrum( time_vector , x_T_LQI , ddx_T_LQI, f_vector , 1);
[picos_ddx_L_LQI , picos_x_L_LQI] = ResponseSpectrum( time_vector , x_L_LQI , ddx_L_LQI, f_vector , 1);

% Computing Response spectra of Adapted
[picos_ddx_T_acq_0  , picos_x_T_acq_0 ] = ResponseSpectrum( time_vector , x_T_acq_0 , ddx_T_acq_0, f_vector , 1);
[picos_ddx_T_acq_1  , picos_x_T_acq_1 ] = ResponseSpectrum( time_vector , x_T_acq_1 , ddx_T_acq_1, f_vector , 1);
[picos_ddx_T_acq_2  , picos_x_T_acq_2 ] = ResponseSpectrum( time_vector , x_T_acq_2 , ddx_T_acq_2, f_vector , 1);
[picos_ddx_L_acq_0  , picos_x_L_acq_0 ] = ResponseSpectrum( time_vector , x_L_acq_0 , ddx_L_acq_0, f_vector , 1);
[picos_ddx_L_acq_1  , picos_x_L_acq_1 ] = ResponseSpectrum( time_vector , x_L_acq_1 , ddx_L_acq_1, f_vector , 1);
[picos_ddx_L_acq_2  , picos_x_L_acq_2 ] = ResponseSpectrum( time_vector , x_L_acq_2 , ddx_L_acq_2, f_vector , 1);


opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;
opts1.PhaseVisible='off'; opts1.YLim={[-30 10]};

fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);ylim([0 ceil(max(picos_ddx_T_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_T_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_T_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',   mean((picos_ddx_tgt_T-picos_ddx_T_LQI).^2 )));
plot(f_vector, picos_ddx_T_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_0).^2 )));
plot(f_vector, picos_ddx_T_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_1).^2 )));
plot(f_vector, picos_ddx_T_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_T_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned).^2 )));%- Normal
plot(f_vector, picos_x_T_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',  mean((picos_x_tgt_T-picos_x_T_LQI).^2 )));
plot(f_vector, picos_x_T_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_0).^2 )));
plot(f_vector, picos_x_T_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_1).^2 )));
plot(f_vector, picos_x_T_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_2).^2 )));

% Create Figures - Longitudinal
fig9 = figure(9);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Parallel');xlim([1 20]);ylim([0 ceil(max(picos_ddx_L_tuned(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Parallel');xlim([0.1 5]);

figure(fig9); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_L,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_L_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF -  MSE= %.2e',      mean((picos_ddx_tgt_L-picos_ddx_L_tuned).^2 )));% - Normal
plot(f_vector, picos_ddx_L_LQI,'--', 'LineWidth' , 2 , 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',   mean((picos_ddx_tgt_L-picos_ddx_L_LQI).^2 )));
plot(f_vector, picos_ddx_L_acq_0 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_0).^2 )));
plot(f_vector, picos_ddx_L_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_1).^2 )));
plot(f_vector, picos_ddx_L_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_L-picos_ddx_L_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_L,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_L_tuned,'--', 'LineWidth' , 2, 'Color', color2, 'DisplayName',sprintf( 'Tuned PIDF - MSE= %.2e',     mean((picos_x_tgt_L-picos_x_L_tuned).^2 )));%- Normal
plot(f_vector, picos_x_L_LQI,'--', 'LineWidth' , 2, 'Color', color3, 'DisplayName',sprintf( 'Optimal Control - MSE= %.2e',  mean((picos_x_tgt_L-picos_x_L_LQI).^2 )));
plot(f_vector, picos_x_L_acq_0, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 0 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_0).^2 )));
plot(f_vector, picos_x_L_acq_1, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 1 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_1).^2 )));
plot(f_vector, picos_x_L_acq_2, '-', 'LineWidth' , 2,  'DisplayName',sprintf( 'Adapted driver 2 - MSE= %.2e', mean((picos_x_tgt_L-picos_x_L_acq_2).^2 )));

%% Simulation using updated driver 3
% LTF_to_TXT_then_load( [ name, '_3.DRV' ] ,'InputFolder',folder)
% x_T_acq_3 = lsim(G_xT_xref ,  x_drv_T_3 , time_vector,'zoh');
% ddx_T_acq_3 = secondDerivativeTime(x_T_acq_3 , t_step);
% % writeTXT_then_LTF(time_vector,x_T_acq_3,ddx_T_acq_3,folder,[ name, '_3.ACQ.txt' ]);
% x_L_acq_3 = lsim(G_xT_xref ,  x_drv_L_3 , time_vector,'zoh');
% ddx_L_acq_3 = secondDerivativeTime(x_L_acq_3 , t_step);
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