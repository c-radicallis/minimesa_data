clear;clc;close all; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-left
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
Ts_fpga= 1/5000;

% spa settings
freq_resolution = 0.1;
win_size = 1/(freq_resolution*Ts); %frequency resolution = 2*pi/(win_size*Ts) [rad/s] = 1/(win_size *Ts)[Hz]
f_vector_full= logspace( log10(2*freq_resolution*2*pi) , log10(100*2*pi) , 100);
f_vector_OL = logspace( log10(2*freq_resolution*2*pi) , log10(100*2*pi) , 20);
f_vector_CL = logspace( log10(2*freq_resolution*2*pi) , log10(100*2*pi) , 20); %linspace( 2*2*freq_resolution*2*pi , 100*2*pi , 6 );

% tfest settings
tfest_opt_OL = tfestOptions('InitialCondition','zero');
np_OL = 6; %number of poles for tf est

tfest_opt_CL = tfestOptions('InitialCondition','zero','EnforceStability',1);
%(30hz,-22.4dB) & (50.1Hz,-37.5dB)   (60.5Hz ,-41.9dB) & (100hz,-71.6dB)
% m_3050 = (-37.5+22.4)/(50-30)*100
% m_50100 = (-71.6+37.5)/(100-50)*100 % *100Hz/dec
% m_30100=(-71.6+22.4)/(100-30)*100
% ≈80 db/dec = 4 polos
np_CL = 4; 

%n4sid settings % SSARX allows unbiased estimates when using closed loop data
nx=5;n4sidOpt = n4sidOptions;n4sidOpt.N4Weight = 'SSARX'; n4sidOpt.Focus = 'simulation';n4sidOpt.InitialState = 'zero';

% ssest Options
ssestOptions = ssestOptions('InitializeMethod','n4sid','EnforceStability',1);

% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[freq_resolution 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;%opts1.Ylim={[-40 10]};

% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert t  o mm

%% Data 12
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
true_tune_12 = pid(7,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data12_OL = iddata(x_acq_T, sv2_acq, Ts);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data12_CL =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);

%% Inputs: u, y, fs, sys  (u,y are column vectors)
nfft    = 2048;
window  = hamming(512);
noverlap= 256;

% data-based coherence
[Cxy_y_u, f_yu] = mscohere(x_acq_T, sv2_acq, window, noverlap, nfft);
[Cxy_y_ref, f_yref] = mscohere(x_acq_T, x_drv_T_0 , window, noverlap, nfft);

figure; hold on;
plot(f_yu, Cxy_y_u, 'LineWidth', 1.5);
plot(f_yref, Cxy_y_ref, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Coherence');
title('Data-based coherence (mscohere)');
grid on;
ylim([0 1]);
legend('OL = x_T / u' , 'CL=x_T / x_{ref}')