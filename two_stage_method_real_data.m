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

%% input file - sine sweep - ddx=1200
sineSweep_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\sineSweep\ddx=1200\';
file = 'sineSweep_ddx=1200_f=1e-5to40.ltf'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', sineSweep_folder , 'OutputFolder', sineSweep_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert to mm

%%  Data sine  - P15
file = 'sineSweep_ddx=1200_f=1e-5to40_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
[~,time_acq_aligned] = alignTimeVectorEnds(time_drv_0, time_acq);



%%  Data sine  - P7
file = 'sineSweep_ddx=1200_f=1e-5to40_P7.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , sineSweep_folder , 'OutputFolder', sineSweep_folder );
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped
[~,time_acq_aligned] = alignTimeVectorEnds(time_drv_0, time_acq);

% step1
S_est = polyest(x_drv_T_0,sv2_acq ,  [[0 15 0 0 0] 0],'Ts',Ts)%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
figure; hold on;
%bodeplot(S0,'r', opts1);
bodeplot(S_est,'g', opts1);
legend; grid on;

% step 2
u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);

figure, hold on;
plot(time_acq_aligned , sv2_acq , 'DisplayName', 'sv2_acq');
plot(time_drv_0,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
legend; grid on;

% step 3
G_est = tfest(u_r_est , x_acq_T , 4 , 'Ts' , Ts)
G_est_to_CL = feedback(15*G_est, 1);
G_direct = tfest( sv2_acq , x_acq_T , 4 , 'Ts' , Ts)
G_direct_to_CL = feedback(15*G_direct, 1);
G_CL = tfest(x_drv_T_0 , x_acq_T , 4 , 'Ts' , Ts)
G_indirect = G_CL/(15*(1-G_CL))

figure;hold on;
bodeplot(G_CL,'y', opts1);
bodeplot(G_est_to_CL , 'g', opts1);
bodeplot( G_direct_to_CL , 'b' , opts1);
bodeplot(G_est , 'g--', opts1);
bodeplot( G_direct , 'b--' , opts1);
bodeplot( G_indirect , 'r--' , opts1);
legend; grid on;

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert t  o mm
clear x_drv_L_0  x_drv_V_0

%%  data_P5
% folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
% file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
% LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
% 
% x_acq_T = x_acq_T*1e3;
% sv2_acq = bits2mm( -sv2_acq ); %output is inverted because the wiring is fliped
% 
% % step1  
% S_est = polyest(x_drv_T_0(1:nmin),sv2_acq ,  [[0 15 0 0 0] 0],'Ts',Ts)%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
% figure; hold on;
% %bodeplot(S0,'r', opts1);
% bodeplot(S_est,'g', opts1);
% legend; grid on;
% 
% % step 2
% u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);
% 
% figure, hold on;
% plot(time_acq , sv2_acq , 'DisplayName', 'sv2_acq');
% plot(time_drv_0,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
% legend; grid on;
% 
% % step 3
% G_est = tfest(u_r_est(1:nmin) , x_acq_T , 4 , 'Ts' , Ts)
% G_est_to_CL = feedback(5*G_est, 1);
% G_direct = tfest( sv2_acq , x_acq_T , 4 , 'Ts' , Ts)
% G_direct_to_CL = feedback(5*G_direct, 1);
% G_CL = tfest(x_drv_T_0(1:nmin) , x_acq_T , 4 , 'Ts' , Ts)
% G_indirect = G_CL/(5*(1-G_CL))
% 
% figure;hold on;
% bodeplot(G_CL,'y', opts1);
% bodeplot(G_est_to_CL , 'g', opts1);
% bodeplot( G_direct_to_CL , 'b' , opts1);
% bodeplot(G_est , 'g--', opts1);
% bodeplot( G_direct , 'b--' , opts1);
% bodeplot( G_indirect , 'r--' , opts1);
% legend; grid on;

%%  data_P7
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm( -sv2_acq ); %output is inverted because the wiring is fliped

% step1
S_est = polyest(x_drv_T_0,sv2_acq ,  [[0 15 0 0 0] 0],'Ts',Ts)%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
figure; hold on;
%bodeplot(S0,'r', opts1);
bodeplot(S_est,'g', opts1);
legend; grid on;

% step 2
u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);

figure, hold on;
plot(time_acq , sv2_acq , 'DisplayName', 'sv2_acq');
plot(time_drv_0,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
legend; grid on;

% step 3
G_est = tfest(u_r_est , x_acq_T , 4 , 'Ts' , Ts)
G_est_to_CL = feedback(7*G_est, 1);
G_direct = tfest( sv2_acq , x_acq_T , 4 , 'Ts' , Ts)
G_direct_to_CL = feedback(7*G_direct, 1);
G_CL = tfest(x_drv_T_0 , x_acq_T , 4 , 'Ts' , Ts)
G_indirect = G_CL/(7*(1-G_CL))

figure;hold on;
bodeplot(G_CL,'y', opts1);
bodeplot(G_est_to_CL , 'g', opts1);
bodeplot( G_direct_to_CL , 'b' , opts1);
bodeplot(G_est , 'g--', opts1);
bodeplot( G_direct , 'b--' , opts1);
bodeplot( G_indirect , 'r--' , opts1);
legend; grid on;

%%  Data P10
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);

x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped

% step1
S_est = polyest(x_drv_T_0(1:nmin),sv2_acq ,  [[0 15 0 0 0] 0],'Ts',Ts)%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
figure; hold on;
%bodeplot(S0,'r', opts1);
bodeplot(S_est,'g', opts1);
legend; grid on;

% step 2
u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);

figure, hold on;
plot(time_acq , sv2_acq , 'DisplayName', 'sv2_acq');
plot(time_drv_0,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
legend; grid on;

% step 3
G_est = tfest(u_r_est(1:nmin) , x_acq_T , 4 , 'Ts' , Ts)
G_est_to_CL = feedback(10*G_est, 1);
G_direct = tfest( sv2_acq , x_acq_T , 4 , 'Ts' , Ts)
G_direct_to_CL = feedback(10*G_direct, 1);
G_CL = tfest(x_drv_T_0(1:nmin) , x_acq_T , 4 , 'Ts' , Ts)
G_indirect = G_CL/(10*(1-G_CL))

figure;hold on;
bodeplot(G_CL,'y', opts1);
bodeplot(G_est_to_CL , 'g', opts1);
bodeplot( G_direct_to_CL , 'b' , opts1);
bodeplot(G_est , 'g--', opts1);
bodeplot( G_direct , 'b--' , opts1);
bodeplot( G_indirect , 'r--' , opts1);
legend; grid on;

%%  Data P15
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);

x_acq_T = x_acq_T*1e3;
sv2_acq = bits2mm(-sv2_acq); %output is inverted because the wiring is fliped

% step1
S_est = polyest(x_drv_T_0,sv2_acq ,  [[0 15 0 0 0] 0],'Ts',Ts)%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
figure; hold on;
%bodeplot(S0,'r', opts1);
bodeplot(S_est,'g', opts1);
legend; grid on;

% step 2
u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);

figure, hold on;
plot(time_acq , sv2_acq , 'DisplayName', 'sv2_acq');
plot(time_drv_0,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
legend; grid on;

% step 3
G_est = tfest(u_r_est , x_acq_T , 4 , 'Ts' , Ts)
G_est_to_CL = feedback(15*G_est, 1);
G_direct = tfest( sv2_acq , x_acq_T , 4 , 'Ts' , Ts)
G_direct_to_CL = feedback(15*G_direct, 1);
G_CL = tfest(x_drv_T_0 , x_acq_T , 4 , 'Ts' , Ts)
G_indirect = G_CL/(15*(1-G_CL))

figure;hold on;
bodeplot(G_CL,'y', opts1);
bodeplot(G_est_to_CL , 'g', opts1);
bodeplot( G_direct_to_CL , 'b' , opts1);
bodeplot(G_est , 'g--', opts1);
bodeplot( G_direct , 'b--' , opts1);
bodeplot( G_indirect , 'r--' , opts1);
legend; grid on;




