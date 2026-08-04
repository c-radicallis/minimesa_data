clear;clc;close all; setappdata(0, 'AutoStagger_LRDown_Last', []);  set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize); 
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
addpath ('C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model' , ...
    'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\' , ...
    'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\functions_matlab')
Ts = 0.005;
Ts_fpga= 1/5000;

%% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[1 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360; 

%% Control channel AI2 Displacement - 16 bit signed integer to mm conversion
a = 0.000485;
b = -0.2;
bits2mm = @(bits) a*bits+b;
mm2bits = @(mm) (mm-b)/a;
clear a b
%%
fir_np=100;
np_CL=4;
np_OL=4;

%% input file - sine sweep - A = 4
% sineSweep_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\sineSweep\A=4_drv\';
% file = 'sineSweep_A=4_f=10e-4to40_CosTaper5percent_0.drv'; 
% LTF_to_TXT_then_load( file , 'InputFolder', sineSweep_folder , 'OutputFolder', sineSweep_folder); % load input drv
% x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
% 
% %  Data sine  - P5
% file = 'sineSweep_A=4_f=10e-4to40_CosTaper5percent_P5.acq';
% LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
% x_acq_T = x_acq_T*1e3;
% sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% Kp=5
% results_P15_sineSweep = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);
% 
% %  Data sine  - P7
% file = 'sineSweep_A=4_f=10e-4to40_CosTaper5percent_P7.acq';
% LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
% x_acq_T = x_acq_T*1e3;
% sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% Kp=7
% results_P7_sineSweep = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

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
results_P7_sineSweep = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%%  Data sine  - P15
file = 'sineSweep_ddx=1200_f=1e-5to40_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=15
results_P15_sineSweep = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
clear x_drv_L_0  x_drv_V_0
 % fs=200; figure; pspectrum(x_drv_T_0,fs); xscale log;

%%  data_P5
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm( -sv2_acq ); %output is inverted because the wiring is fliped
Kp=5
results_P5_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%%  data_P7
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm( -sv2_acq ); %output is inverted because the wiring is fliped
Kp=7
results_P7_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%%  Data P10
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=10
results_P10_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
Kp=15
results_P15_pink = twoStageMethod(Kp , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%%  Data wc_10Hz
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\pink_noise';
file = 'wc_10Hz_0.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201 , 'OutputFolder', folder_1201);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% w_c = 10Hz
Kp = 6.75;
Ki = 68;
Kd = .0712891;
Tf = 0.00161743;
controller_at_200 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
%controller_at_200 = pid(Kp,Ki,Kd,Tf,Ts)
results_wc10_pink = twoStageMethod( controller_at_200 , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);

%%  Data wc_15Hz
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\pink_noise';
file = 'wc_15Hz_0.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201 , 'OutputFolder', folder_1201);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% w_c = 15Hz
Kp = 8.75
Ki = 139
Kd = 0.125977
Tf = 0.00109863
controller_at_200 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
%controller_at_200 = pid(Kp,Ki,Kd,Tf,Ts)
results_wc15_pink = twoStageMethod(controller_at_200 , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);


%%  Data wc_20Hz
folder_1201 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\12-1-2026\pink_noise';
file = 'wc_20Hz_0.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_1201 , 'OutputFolder', folder_1201);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% w_c = 20Hz
Kp = 8.5625;
Ki = 75;
Kd = 0.18457;
Tf = 0.000167847;
controller_at_200 = d2d(pid(Kp,Ki,Kd,Tf,Ts_fpga),Ts)
%controller_at_200 = pid(Kp,Ki,Kd,Tf,Ts)
results_wc20_pink = twoStageMethod(controller_at_200 , fir_np, np_CL , np_OL,  Ts , opts1, sv2_acq, x_drv_T_0, time_drv_0, time_acq, x_acq_T,1);


%% Let's compare open loop tranfers fucntions
%close all;
opts1.PhaseVisible='off'; opts1.XLim={[0.1 40]};  opts1.YLim={[-40 5]};
full_legend = {'P7 (sweep)', 'P15 (sweep)',  'P5','P7', 'P10', 'P15', 'wc10', 'wc15', 'wc20'};
nice_legend = {'P7', 'P10', 'P15',  'wc15', 'wc20'};
exc_sweep_legend =  {  'P5','P7', 'P10', 'P15', 'wc10', 'wc15', 'wc20'};

figure; hold on; opts1.Title.String='OL_{direct}';
bodeplot(results_P7_sineSweep.OL_direct, opts1);
bodeplot(results_P15_sineSweep.OL_direct, opts1);
bodeplot(results_P5_pink.OL_direct, opts1);
bodeplot(results_P7_pink.OL_direct, opts1);
bodeplot(results_P10_pink.OL_direct, opts1);
bodeplot(results_P15_pink.OL_direct, opts1);
bodeplot(results_wc10_pink.OL_direct, opts1);
bodeplot(results_wc15_pink.OL_direct, opts1);
bodeplot(results_wc20_pink.OL_direct, opts1);
legend(full_legend,"Location",'southwest'); grid on;

figure; hold on; opts1.Title.String='OL_{indirect}';
bodeplot(results_P7_sineSweep.OL_indirect, opts1);
bodeplot(results_P15_sineSweep.OL_indirect, opts1);
bodeplot(results_P5_pink.OL_indirect, opts1);
bodeplot(results_P7_pink.OL_indirect, opts1);
bodeplot(results_P10_pink.OL_indirect, opts1);
bodeplot(results_P15_pink.OL_indirect, opts1);
bodeplot(results_wc10_pink.OL_indirect, opts1);
bodeplot(results_wc15_pink.OL_indirect, opts1);
bodeplot(results_wc20_pink.OL_indirect, opts1);
legend(full_legend,"Location",'southwest'); grid on;

figure; hold on; opts1.Title.String='OL_{two-stage}';
bodeplot(results_P7_sineSweep.OL_two_stage, opts1);
bodeplot(results_P15_sineSweep.OL_two_stage, opts1);
bodeplot(results_P5_pink.OL_two_stage, opts1);
bodeplot(results_P7_pink.OL_two_stage, opts1);
bodeplot(results_P10_pink.OL_two_stage, opts1);
bodeplot(results_P15_pink.OL_two_stage, opts1);
bodeplot(results_wc10_pink.OL_two_stage, opts1);
bodeplot(results_wc15_pink.OL_two_stage, opts1);
bodeplot(results_wc20_pink.OL_two_stage, opts1);
legend(full_legend,"Location",'southwest'); grid on;

% Let's compare CLOSED LOOP tranfers fucntions
opts1.XLim={[1 40]}; opts1.YLim={[-25 8]};

% figure; hold on; opts1.Title.String='CL from OL_{direct}';
% % % bodeplot(results_P7_sineSweep.CL_from_OL_direct, opts1);
% % % bodeplot(results_P15_sineSweep.CL_from_OL_direct, opts1);
% % bodeplot(results_P5_pink.CL_from_OL_direct, opts1);
% bodeplot(results_P7_pink.CL_from_OL_direct, opts1);
% bodeplot(results_P10_pink.CL_from_OL_direct, opts1);
% bodeplot(results_P15_pink.CL_from_OL_direct, opts1);
% % % bodeplot(results_wc10_pink.CL_from_OL_direct, opts1);
% % bodeplot(results_wc15_pink.CL_from_OL_direct, opts1);
% % bodeplot(results_wc20_pink.CL_from_OL_direct, opts1);
% % legend('P7','P10','P15','wc10','wc15','wc20'); grid on;

% figure; hold on; opts1.Title.String='CL';
% % bodeplot(results_P7_sineSweep.CL, opts1);
% % bodeplot(results_P15_sineSweep.CL, opts1);
% bodeplot(results_P5_pink.CL, opts1);
% bodeplot(results_P7_pink.CL, opts1);
% bodeplot(results_P10_pink.CL, opts1);
% bodeplot(results_P15_pink.CL, opts1);
% bodeplot(results_wc10_pink.CL, opts1);
% bodeplot(results_wc15_pink.CL, opts1);
% bodeplot(results_wc20_pink.CL, opts1);
% legend(exc_sweep_legend,"Location",'southwest'); grid on;

% figure; hold on; opts1.Title.String='CL from OL_{est nonLin}';
% % % bodeplot(results_P7_sineSweep.CL_from_OL_two_stage, opts1);
% % % bodeplot(results_P15_sineSweep.CL_from_OL_two_stage, opts1);
% % bodeplot(results_P5_pink.CL_from_OL_two_stage, opts1);
% bodeplot(results_P7_pink.CL_from_OL_two_stage, opts1);
% bodeplot(results_P10_pink.CL_from_OL_two_stage, opts1);
% bodeplot(results_P15_pink.CL_from_OL_two_stage, opts1);
% % % bodeplot(results_wc10_pink.CL_from_OL_two_stage, opts1);
% % bodeplot(results_wc15_pink.CL_from_OL_two_stage, opts1);
% % bodeplot(results_wc20_pink.CL_from_OL_two_stage, opts1);
% % legend('P7','P10','P15','wc10','wc15','wc20'); grid on; %
% 

%%
% FolderName = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\id_results';   % Your destination folder
% baseFolder = FolderName;   % Base folder where you want to create the timestamped subfolder
% ts = datestr(now, 'yyyymmdd_HHMM');  % Create a timestamp string, e.g. '20250709_1530'
% timeDir = fullfile(baseFolder, ts);  % Build the full path to the new folder
% if ~exist(timeDir, 'dir')% Create it if it doesn't already exist
%     mkdir(timeDir)
% end
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
%   savefig(fullfile(timeDir, [FigName '.fig']));
%   exportgraphics(FigHandle, fullfile(timeDir, [FigName '.jpg']));
% end
