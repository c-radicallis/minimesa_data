addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);

Ts = 0.005;
Ts_fpga= 1/5000;

% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[freq_resolution 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;%opts1.Ylim={[-40 10]};

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert t  o mm
clear x_drv_L_0  x_drv_V_0

%%  data_P7
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
nmin = min(numel(time_drv_0), numel(time_acq));

x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped

erro = x_drv_T_0(1:nmin)-x_acq_T(1:nmin);
data_P7_r_to_u = iddata( sv2_acq(1:nmin) , erro , Ts ); data_P7_r_to_u.InputName  = 'erro';data_P7_r_to_u.OutputName = 'sv2_acq';data_P7_r_to_u.TimeUnit   = 'seconds';
% figure; plot(data_P7_r_to_u);
