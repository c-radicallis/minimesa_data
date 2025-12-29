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

fir_np=40;
np_CL=4;
np_OL=4;

%% input file - sine sweep - ddx=1200
sineSweep_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\sineSweep\ddx=1200\';
file = 'sineSweep_ddx=1200_f=1e-5to40.ltf'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', sineSweep_folder , 'OutputFolder', sineSweep_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm


%%  Data sine  - P7
file = 'sineSweep_ddx=1200_f=1e-5to40_P7.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=7
results_P7_sineSweep = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%%  Data sine  - P15
file = 'sineSweep_ddx=1200_f=1e-5to40_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=15
results_P15_sineSweep = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
clear x_drv_L_0  x_drv_V_0

%%  data_P5
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm( -sv2_acq ); %output is inverted because the wiring is fliped
Kp=5
results_P5_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%%  data_P7
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm( -sv2_acq ); %output is inverted because the wiring is fliped
Kp=7
results_P7_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);


%%  Data P10
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=10
results_P10_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);

%%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=15
results_P15_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T);