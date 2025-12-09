clear;clc;close all; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-left
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
Ts_fpga= 1/5000;

% % spa settings
freq_resolution = 0.1;
win_size = 1/(freq_resolution*Ts); %frequency resolution = 2*pi/(win_size*Ts) [rad/s] = 1/(win_size *Ts)[Hz]

% tfest settings
tfest_opt_OL = tfestOptions('InitialCondition','auto');
np_OL = 4; %number of poles for tf est

tfest_opt_CL = tfestOptions('InitialCondition','auto','EnforceStability',1);
np_CL = 4; 

np_C_est = 20;

%n4sid settings % SSARX allows unbiased estimates when using closed loop data
% nx=5;n4sidOpt = n4sidOptions;n4sidOpt.N4Weight = 'SSARX'; n4sidOpt.Focus = 'simulation';n4sidOpt.InitialState = 'zero';

% ssest Options
% ssestOptions = ssestOptions('InitializeMethod','n4sid','EnforceStability',1);

% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[freq_resolution 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;%opts1.Ylim={[-40 10]};

% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert t  o mm
clear x_drv_L_0  x_drv_V_0

% Control channel AI2 Displacement - 16 bit signed integer to mm conversion
a = 0.000485;
b = -0.2;
bits2mm = @(bits) a*bits+b;
mm2bits = @(mm) (mm-b)/a;
clear a b

fcut=30; % Hz - cutoff frequency for model order reduction

% %%  data_P5
% folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
% file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
% nmin = min(numel(time_drv_0), numel(time_acq));
% 
% x_acq_T = x_acq_T*1e3;
% sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped % converted to mm
% 
% erro = x_drv_T_0(1:nmin)-x_acq_T(1:nmin);
% data_P5_r_to_u = iddata( sv2_acq(1:nmin) , erro , Ts ); data_P5_r_to_u.InputName  = 'erro';data_P5_r_to_u.OutputName = 'sv2_acq';data_P5_r_to_u.TimeUnit   = 'seconds';
% % figure; plot(data_P5_r_to_u);
% 
% C_estimated = tfest(data_P5_r_to_u ,40 ,'Ts',Ts,tfest_opt_OL)
% %[C_estimated_reduced , ~]=freqsep( C_estimated , deg2rad(fcut));
% u_simulated = lsim(C_estimated , erro , time_acq);
% figure;hold on;title('Controler:  K=5');title('Controler:  K=5');
% plot(time_acq , sv2_acq,'DisplayName' , 'sv2_acq');
% plot(time_acq, u_simulated); %, 'DisplayName', ['static gain K = ' num2str(C_estimated.Numerator,'%.3f')]);
% grid on; legend;

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

C_estimated = tfest(data_P7_r_to_u ,np_C_est  ,'Ts',Ts,tfest_opt_OL)
% C_estimated_nlarx = nlarx(data_P7_r_to_u,[10 10 1])
figure;compare( data_P7_r_to_u , C_estimated ); %, C_estimated_nlarx);
% C_estimated_1p = tfest(data_P7_r_to_u , 1 ,'Ts',Ts,tfest_opt_OL)
% u_simulated = lsim(C_estimated , erro , time_acq);
% u_simulated_nlarx = lsim( C_estimated_nlarx , erro , time_acq);
% figure;hold on;title('Controler:  K=7');
% plot(time_acq , sv2_acq,'DisplayName' , 'sv2_acq');
% plot(time_acq, u_simulated); %, 'DisplayName', ['static gain K = ' num2str(C_estimated.Numerator,'%.3f')]);
% plot(time_acq , u_simulated_nlarx,'DisplayName' , 'nlarx');
% grid on;legend

% %%  Data P10
% folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
% file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
% nmin = min(numel(time_drv_0), numel(time_acq));
% 
% x_acq_T = x_acq_T*1e3;
% sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
% 
% erro = x_drv_T_0(1:nmin)-x_acq_T(1:nmin);
% data_P10_r_to_u = iddata( sv2_acq(1:nmin) , erro , Ts ); data_P10_r_to_u.InputName  = 'erro';data_P10_r_to_u.OutputName = 'sv2_acq';data_P10_r_to_u.TimeUnit   = 'seconds';
% % figure; plot(data_P10_r_to_u);
% 
% C_estimated = tfest(data_P10_r_to_u ,40 ,'Ts',Ts,tfest_opt_OL)
% % C_estimated_1p = tfest(data_P10_r_to_u , 1 ,'Ts',Ts,tfest_opt_OL)
% u_simulated = lsim(C_estimated , erro , time_acq);
% % u_simulated_1p = lsim(% C_estimated_1p , erro , time_acq);
% figure;hold on;title('Controler:  K=10');
% plot(time_acq , sv2_acq,'DisplayName' , 'sv2_acq');
% plot(time_acq, u_simulated); %, 'DisplayName', ['static gain K = ' num2str(C_estimated.Numerator,'%.3f')]);
% %plot(time_acq , u_simulated_1p,'DisplayName' , '1 pole');
% grid on;legend

%%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
nmin = min(numel(time_drv_0), numel(time_acq));

x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped

erro = x_drv_T_0(1:nmin)-x_acq_T(1:nmin);
data_P15_r_to_u = iddata( sv2_acq(1:nmin) , erro , Ts ); data_P15_r_to_u.InputName  = 'erro';data_P15_r_to_u.OutputName = 'sv2_acq';data_P15_r_to_u.TimeUnit   = 'seconds';
% figure; plot(data_P15_r_to_u);

C_estimated = tfest(data_P15_r_to_u ,np_C_est ,'Ts',Ts,tfest_opt_OL)
% C_estimated_nlarx = nlarx(data_P15_r_to_u,[10 10 1])
figure;compare( data_P15_r_to_u , C_estimated ); %, C_estimated_nlarx);
% % C_estimated_1p = tfest(data_P15_r_to_u , 1 ,'Ts',Ts,tfest_opt_OL)
% u_simulated = lsim(C_estimated , erro , time_acq);
% % u_simulated_1p = lsim(% C_estimated_1p , erro , time_acq);
% figure;hold on;title('Controler:  K=15');
% plot(time_acq , sv2_acq,'DisplayName' , 'sv2_acq');
% plot(time_acq, u_simulated); %, 'DisplayName', ['static gain K = ' num2str(C_estimated.Numerator,'%.3f')]);
% %plot(time_acq , u_simulated_1p,'DisplayName' , '1 pole');
% grid on;legend
