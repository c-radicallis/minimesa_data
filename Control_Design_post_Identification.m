% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data2 in MATLAB with the System Identification Toolbox.

% clear;clc;
% addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
% func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
% addpath(func_folder);
% Ts = 0.005;
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[6 50]}; opts1.PhaseMatching='on'; opts1.Grid='on';
clc;
close all;

OL_direct = results_P15_pink.OL_direct
OL_indirect = results_P15_pink.OL_indirect
OL_est_nonLin = results_P15_pink.OL_est_nonLin
CL_from_OL_direct = results_P15_pink.CL_from_OL_direct
CL = results_P15_pink.CL
CL_from_OL_est_nonLin = results_P15_pink.CL_from_OL_est_nonLin

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

%%
G_open=OL_est_nonLin;
G_closed = results_P15_pink.CL_from_OL_est_nonLin;

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

%%
% obs = vpa(obsv(G_open));
% r_obsv = rank(obs)
% ctrlb = vpa(ctrb(G_open));
% r_ctrlb = rank(ctrlb)

G_open = ss(OL_est_nonLin);

n_states=size(G_open.A,1);
G_open.StateName  = arrayfun(@(k) sprintf('x%d',k), 1:n_states, 'UniformOutput', false);

plant_aug = ss(G_open.A, G_open.B,[eye(n_states);G_open.C],0 , Ts);%_fpga);
plant_aug.InputName = {'i_sv'};   % plant input: control signal
plant_aug.OutputName = [G_open.StateName ; {'y_xT'}];  % plant output

sumblk1 = sumblk('e = x_ref - y_xT'); % Compute the error signal: e = r - y

integrator = tf(1,[1 0], Ts);%_fpga); % The integrator integrates the tracking error.
integrator.InputName = {'e'};    % error: e = r - y
integrator.OutputName = {'xi'};  % integrated error

Q = diag([0*ones(1,n_states),1e5]);%1e3*diag([zeros(1,n_states),1]);%blkdiag(eye(nx), eye(ny));
R = eye(size(G_open.B,2));
K_lqi = lqi(G_open, Q, R)% Design the LQI controller for the original system

K  = K_lqi(1:n_states);      % state feedback gains
Ki = K_lqi(end);        % integrator gain
controller = ss([], [], [], -[K, Ki]); %   u = -[K  Ki] * [x; xi]
controller.InputName = [ G_open.StateName ; {'xi'}];
controller.OutputName = {'i_sv'};
Optimal_closed_loop = minreal(connect(plant_aug,  controller , integrator, sumblk1, 'x_ref','y_xT'))

%figure;pzmap(Optimal_closed_loop);

%
fig9 = figure(9);ax9 = axes(fig9); hold(ax9, 'on');
bodeplot(G_closed,opts1);
bodeplot(G_PIDF_10Hz,opts1);    
bodeplot(G_PIDF_15Hz,opts1); 
bodeplot(G_PIDF_20Hz,opts1);    
bodeplot(Optimal_closed_loop,opts1);
legend();


%% Finding Response Spectre of Ground
% f_i=0.1; %freq inicial
% f_n=30;  %freq final
% n_points = 2e2;
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
% x_Optimal = lsim(d2d(Optimal_closed_loop,Ts) ,  x_ref , t_vector,'zoh');
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