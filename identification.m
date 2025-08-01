% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data in MATLAB with the System Identification Toolbox.
clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);

% Load input drv
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
addpath(input_file_folder);
file = 'pink_noise_40Hz_T3mm_0.drv';
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);


% load output acq
file = 'pink_noise_40Hz_T3mm_scl=1_0.acq';
scale = 1;

x_drv_T_0=scale*x_drv_T_0;
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);
Ts = time_acq(2);
valdata = iddata(x_acq_T, x_drv_T_0, Ts);
valdata.InputName  = 'x_drv_T_0';
valdata.OutputName = 'x_acq_T';
valdata.TimeUnit   = 'seconds';


% load output acq
file = 'pink_noise_40Hz_T3mm_scl=1.2_0.acq';
scale = 1.2;

x_drv_T_0=scale*x_drv_T_0;
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder);
int_int_ddx_acq_T = lsim(tf(1,[1 0 0]) , ddx_acq_T,time_acq);

%%
% --- Step 2: Create an iddata object ---
% Determine lengths
n1 = numel(x_drv_T_0);
n2 = numel(x_acq_T);
nmin = min(n1, n2);
% Truncate both to the shortest length
x_drv_T_0 = x_drv_T_0(1:nmin);
x_acq_T = x_acq_T(1:nmin);
ddx_acq_T = ddx_acq_T(1:nmin);

data = iddata(x_acq_T, x_drv_T_0, Ts);
data.InputName  = 'x_drv_T_0';
data.OutputName = 'x_acq_T';
data.TimeUnit   = 'seconds';

% data2 = iddata(int_int_ddx_acq_T, x_drv_T_0, Ts);
% data2.InputName  = 'x_drv_T_0';
% data2.OutputName = 'int_int_ddx_acq_T';
% data2.TimeUnit   = 'seconds';

nz  = 1;
np = 2;
tf_model1 = tfest(data,np,nz)%,'Ts',Ts )

na = 2; nb = 2; nk = 1;  % orders and input delay
arx_model = arx(data, [na nb nk]) % Alternative ARX estimation:

%%
% % Suppose `data` is your iddata for estimation, and `valdata` is separate validation data.
% na = 1:10;         % possible orders for A(z^-1)
% nb = 1:10;         % possible orders for B(z^-1)
% nk = 1:5;          % possible input delay(s)
% 
% % 1) Build the grid of candidate [na nb nk]
% M = struc(na, nb, nk);
% 
% % 2) Compute loss metrics for each candidate
% %    V.Fit contains the chosen loss (default: AIC),
% %    V.FPE the Final Prediction Error, etc.
% V = arxstruc(data, valdata, M);
% 
% [~, bestIdxFPE] = selstruc(V,'mse');
% bestOrders = M(bestIdxFPE,:);
% 
% disp(['Best ARX structure: na=' num2str(bestOrders(1)) ...
%       ', nb=' num2str(bestOrders(2)) ...
%       ', nk=' num2str(bestOrders(3)) ]);
%%

nx = 3;
n4sid_sys = n4sid( x_drv_T_0,x_acq_T,nx,'Ts',Ts)
% Copy names from data into your model:
n4sid_sys.InputName  = data.InputName;
n4sid_sys.OutputName = data.OutputName;


%%
% optCVA = n4sidOptions('N4weight','CVA');
% optSSARX = n4sidOptions('N4weight','SSARX');
% 
% sysCVA = n4sid(data,nx,optCVA);
% sysCVA.InputName  = data.InputName;
% sysCVA.OutputName = data.OutputName;
% sysSSARX = n4sid(data,nx,optSSARX);
% sysSSARX.InputName  = data.InputName;
% sysSSARX.OutputName = data.OutputName;
% 
% figure(3)
% compare(data,sysCVA,sysSSARX);

%%
figure(2);
compare(data,tf_model1,arx_model,n4sid_sys)

figure(3);
compare(valdata,tf_model1,arx_model,n4sid_sys)

% y1 = lsim(tf_model1, x_drv_T_0, time_acq);
% y2 = lsim(arx_model, x_drv_T_0, time_acq);
%plot(time_acq, x_acq_T, 'b', time_acq, y1, 'r--',time_acq, y2, 'm--');
%legend('Measured y','TFest' , 'Arx');
% title('Model Fit');
% xlabel('Time (s)');
% ylabel('Output');

% subplot(3,1,2);
% resid(data, tf_model1);
% title('Residual Analysis');
% subplot(3,1,3);
% resid(data, arx_model);
% title('Residual Analysis');


%%
fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on');
opts1=bodeoptions('cstprefs');
opts1.FreqUnits = 'Hz';
opts1.XLim={[1 40]};
bodeplot(tf_model1,opts1)
bodeplot(arx_model,opts1)
bodeplot(n4sid_sys,opts1)
legend()
grid on


