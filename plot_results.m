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

%%  Tolmezzo tgt
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\';
file = 'TolmezzoReducedScale.tgt'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', folder_1201); % load input drv
scale = 0.16;
x_tgt_T = scale*x_tgt_T; % convert to mm
ddx_tgt_T = scale*ddx_tgt_T; % convert to mm

% Tolmezzo tune 10hz, 15hz , 20hz results
folder_1201_tolmezzo ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\tolmezzo_scl0.16\';
file = 'tune_cutoff_10Hz.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
x_acq_T_10hz = x_acq_T;
%sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
ddx_acq_T_10hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

file = 'tune_cutoff_15hz.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
x_acq_T_15hz= x_acq_T;
%sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
ddx_acq_T_15hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

file = 'tune_cutoff_20hz.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
x_acq_T_20hz = x_acq_T;
%sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
ddx_acq_T_20hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

%%Computing Response spectra of Adapted
folder_1201_tolmezzo_drv_acq ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\tolmezzo_drv&acq\';
file = 'TolmezzoReducedScale_4.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo_drv_acq );
ddx_acq_T_comp = secondDerivativeTime(x_acq_T , Ts);

%% Create Figures - Transversal
close all;

% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);

% Finding Response Spectre  of tuned pidf 10h 15hz 20hz
[picos_ddx_T_tuned_10hz, picos_x_T_tuned_10hz] = ResponseSpectrum( time_vector , x_acq_T_10hz , ddx_acq_T_10hz_comp, f_vector , 1);
[picos_ddx_T_tuned_15hz, picos_x_T_tuned_15hz] = ResponseSpectrum( time_vector , x_acq_T_15hz , ddx_acq_T_15hz_comp, f_vector , 1);
[picos_ddx_T_tuned_20hz, picos_x_T_tuned_20hz] = ResponseSpectrum( time_vector , x_acq_T_20hz , ddx_acq_T_20hz_comp, f_vector , 1);

% Adapted driver response spectra
[picos_ddx_T_acq_4  , picos_x_T_acq_4 ] = ResponseSpectrum( time_vector , x_acq_T, ddx_acq_T_comp, f_vector , 1);

baseFolder = folder_1201;   % Base folder where you want to create the timestamped subfolder
ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
    mkdir(timeDir)
end

fig8 = figure(8);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);%ylim([0 ceil(max(picos_ddx_T_tuned_20hz(1:385,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig8); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_T_tuned_10hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=10Hz -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned_10hz).^2 )));% - Normal
plot(f_vector, picos_ddx_T_tuned_15hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=15Hz -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned_15hz).^2 )));% - Normal
plot(f_vector, picos_ddx_T_tuned_20hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=20Hz -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned_20hz).^2 )));% - Normal
plot(f_vector, picos_ddx_T_acq_4 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 4 - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_4).^2 )));
% plot(f_vector, picos_ddx_T_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_1).^2 )));
% plot(f_vector, picos_ddx_T_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_T_tuned_10hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=10Hz - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned_10hz).^2 )));%- Normal
plot(f_vector, picos_x_T_tuned_15hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=15Hz - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned_15hz).^2 )));%- Normal
plot(f_vector, picos_x_T_tuned_20hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=20Hz - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned_20hz).^2 )));%- Normal
plot(f_vector, picos_x_T_acq_4, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 4 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_4).^2 )));

%%  laquila tgt
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\';
file = 'laquilaReducedScale.tgt'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', folder_1201); % load input drv
scale = 0.18;
x_tgt_T = scale*x_tgt_T; % convert to mm
ddx_tgt_T = scale*ddx_tgt_T; % convert to mm

% laquila tune 10hz, 15hz , 20hz results
folder_1201_tolmezzo ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\laquila_scl0.18\';
file = 'tune_cutoff_10Hz.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
x_acq_T_10hz = x_acq_T;
%sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
ddx_acq_T_10hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

file = 'tune_cutoff_15hz.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
x_acq_T_15hz= x_acq_T;
%sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
ddx_acq_T_15hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

file = 'tune_cutoff_20hz.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo );
x_acq_T_20hz = x_acq_T;
%sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
ddx_acq_T_20hz_comp = secondDerivativeTime(x_acq_T , Ts); % No acceleration data was colected so i'll compute it

%%Computing Response spectra of Adapted
folder_1201_laquila_drv_acq ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\laquila_drv&acq\';
file = 'laquilaReducedScale_4.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201_tolmezzo_drv_acq );
ddx_acq_T_comp = secondDerivativeTime(x_acq_T , Ts);

%% Create Figures 
close all;

% Response Spectra settings
f_i=0.1; %freq inicial
f_n=30;  %freq final
n_points = 5e2;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

% Finding Response Spectre  of Target
[picos_ddx_tgt_T , picos_x_tgt_T] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 1);

% Finding Response Spectre  of tuned pidf 10h 15hz 20hz
[picos_ddx_T_tuned_10hz, picos_x_T_tuned_10hz] = ResponseSpectrum( time_vector , x_acq_T_10hz , ddx_acq_T_10hz_comp, f_vector , 1);
[picos_ddx_T_tuned_15hz, picos_x_T_tuned_15hz] = ResponseSpectrum( time_vector , x_acq_T_15hz , ddx_acq_T_15hz_comp, f_vector , 1);
[picos_ddx_T_tuned_20hz, picos_x_T_tuned_20hz] = ResponseSpectrum( time_vector , x_acq_T_20hz , ddx_acq_T_20hz_comp, f_vector , 1);

% Adapted driver response spectra
[picos_ddx_T_acq_4  , picos_x_T_acq_4 ] = ResponseSpectrum( time_vector , x_acq_T, ddx_acq_T_comp, f_vector , 1);

baseFolder = folder_1201;   % Base folder where you want to create the timestamped subfolder
ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
    mkdir(timeDir)
end

fig9 = figure(9);subplot(121); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 20]);%ylim([0 ceil(max(picos_ddx_T_tuned_20hz(1:395,1))) ])
subplot(122);grid on;xlabel('Frequency (Hz)');ylabel('Displacement (m)');title('Displacement Response Spectra - Fault Normal');xlim([0.1 5]);
color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';% Define colors for lines 1/3 and 2/4

figure(fig9); subplot(121); grid on; legend(); hold on;
plot(f_vector, picos_ddx_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_T_tuned_10hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=10Hz -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned_10hz).^2 )));% - Normal
plot(f_vector, picos_ddx_T_tuned_15hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=15Hz -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned_15hz).^2 )));% - Normal
plot(f_vector, picos_ddx_T_tuned_20hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=20Hz -  MSE= %.2e',      mean((picos_ddx_tgt_T-picos_ddx_T_tuned_20hz).^2 )));% - Normal
plot(f_vector, picos_ddx_T_acq_4 ,'-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 4 - MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_4).^2 )));
% plot(f_vector, picos_ddx_T_acq_1 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 1 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_1).^2 )));
% plot(f_vector, picos_ddx_T_acq_2 ,'-', 'LineWidth' , 2, 'DisplayName',sprintf( 'Adapted driver 2 -  MSE= %.2e', mean((picos_ddx_tgt_T-picos_ddx_T_acq_2).^2 )));

subplot(122); grid on;legend();hold on;
plot(f_vector, picos_x_tgt_T,'-', 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target ');%- Normal
plot(f_vector, picos_x_T_tuned_10hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=10Hz - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned_10hz).^2 )));%- Normal
plot(f_vector, picos_x_T_tuned_15hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=15Hz - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned_15hz).^2 )));%- Normal
plot(f_vector, picos_x_T_tuned_20hz,'--', 'LineWidth' , 2, 'DisplayName',sprintf( 'PIDF ω_c=20Hz - MSE= %.2e',     mean((picos_x_tgt_T-picos_x_T_tuned_20hz).^2 )));%- Normal
plot(f_vector, picos_x_T_acq_4, '-', 'LineWidth' , 2, 'Color', color4, 'DisplayName',sprintf( 'Adapted driver 4 - MSE= %.2e', mean((picos_x_tgt_T-picos_x_T_acq_4).^2 )));



%%
set(fig8, 'WindowState', 'maximized');
exportgraphics(fig8,fullfile(timeDir,'Response_Spectra_N.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');

set(fig9, 'WindowState', 'maximized');
exportgraphics(fig9,fullfile(timeDir,'Response_Spectra_P.png'),'Resolution', 300,'BackgroundColor', 'white','ContentType', 'image');