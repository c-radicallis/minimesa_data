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

%%  Data 11
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
true_tune_11 = pid(5,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data11_OL = iddata(x_acq_T, sv2_acq, Ts);data11_OL.InputName  = 'sv2_acq';data11_OL.OutputName = 'x_acq_T';data11_OL.TimeUnit   = 'seconds';
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data11_CL =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);data11_CL.InputName  = 'x_drv_T_0';data11_CL.OutputName = 'x_acq_T';data11_CL.TimeUnit   = 'seconds';
%
spa_data11_OL_full = spa(data11_OL, win_size, f_vector_full);
spa_data11_OL = spa(data11_OL, win_size, f_vector_OL);
spa_data11_CL_full = spa(data11_CL, win_size, f_vector_full);
spa_data11_CL = spa(data11_CL, win_size, f_vector_CL);

tfest_spa_data11_OL = tfest(spa_data11_OL,np_OL,'Ts',Ts,tfest_opt_OL)
tfest_spa_data11_CL = tfest(spa_data11_CL,np_CL,'Ts',Ts,tfest_opt_CL)
ssest_data11_CL =  ssest(spa_data11_CL_full,6,ssestOptions)
%
spa_OL_from_Tune_and_CL_spa =  spa_data11_CL_full/(d2d(true_tune_11,Ts)*(1-spa_data11_CL_full));%minreal() 
tfest_spa_OL_from_Tune_and_CL_spa = tfest(spa_OL_from_Tune_and_CL_spa,np_OL,'Ts',Ts,tfest_opt_OL)
spa_CL_from_Tune_and_OL_spa = feedback(d2d(true_tune_11,Ts)*spa_data11_OL_full, 1);
CL_from_Tune_and_OL_tfest = feedback(true_tune_11*d2d(tfest_spa_data11_OL,Ts_fpga), 1);

% Open Loop
fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
bodeplot(spa_data11_OL_full,"k.");
bodeplot(spa_data11_OL   ,opts1,"r*");%showConfidence(h)
bodeplot(tfest_spa_data11_OL   ,opts1,"b");%showConfidence(h);
bodeplot(spa_OL_from_Tune_and_CL_spa   ,opts1,"g*");
bodeplot(tfest_spa_OL_from_Tune_and_CL_spa   ,opts1,"g-");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","OL from tune and CL"); grid on;

% Closed Loop
fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop'); 
bodeplot(spa_data11_CL_full,"k.");
bodeplot(spa_data11_CL   ,opts1,"r*");% showConfidence(h)
bodeplot(tfest_spa_data11_CL   ,opts1,"b");% showConfidence(h);
bodeplot(spa_CL_from_Tune_and_OL_spa   ,opts1,"g*"); 
bodeplot(CL_from_Tune_and_OL_tfest   ,opts1,"g-");
bodeplot(ssest_data11_CL,opts1,"y-");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","CL from tune and OL spa","CL from tune and OL tfest"); grid on;

%
Ymodel_tfest_spaOL = lsim(tfest_spa_data11_OL,sv2_acq,time_acq); 
E = Ymodel_tfest_spaOL - x_acq_T;
half=floor(length(time_acq)/2);
figure(91); autocorr(E(1:half),NumLags=300);
% R_XE = xcorr(E,x_acq_T,'coeff'); %  max(E) = 8e+277
% figure(92); plot( -time_acq(end):Ts:time_acq(end) , R_XE)

Ymodel_tfest_spaOL = lsim(tfest_spa_OL_from_Tune_and_CL_spa,sv2_acq,time_acq); 
E = Ymodel_tfest_spaOL - x_acq_T;
half=floor(length(time_acq)/2);
figure(92); autocorr(E(1:half),NumLags=100);

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
%
spa_data12_OL_full = spa(data12_OL, win_size, f_vector_full);
spa_data12_OL = spa(data12_OL, win_size, f_vector_OL);
spa_data12_CL_full = spa(data12_CL, win_size, f_vector_full);
spa_data12_CL = spa(data12_CL, win_size, f_vector_CL);

tfest_spa_data12_OL = tfest(spa_data12_OL,np_OL,'Ts',Ts,tfest_opt_OL)
tfest_spa_data12_CL = tfest(spa_data12_CL,np_CL,'Ts',Ts,tfest_opt_CL)
ssest_data12_CL =  ssest(spa_data12_CL_full,6,ssestOptions)
%
spa_OL_from_Tune_and_CL_spa =  spa_data12_CL_full/(d2d(true_tune_12,Ts)*(1-spa_data12_CL_full));%minreal() 
tfest_spa_OL_from_Tune_and_CL_spa = tfest(spa_OL_from_Tune_and_CL_spa,np_OL,'Ts',Ts,tfest_opt_OL)
CL_from_tune_and_tfest_spa_OL_from_Tune_and_CL_spa = feedback(true_tune_12*d2d(tfest_spa_OL_from_Tune_and_CL_spa,Ts_fpga), 1);

spa_CL_from_Tune_and_OL_spa = feedback(d2d(true_tune_12,Ts)*spa_data12_OL_full, 1);
CL_from_Tune_and_OL_tfest = feedback(true_tune_12*d2d(tfest_spa_data12_OL,Ts_fpga), 1);

% Open Loop
fig12 = figure(12);ax12 = axes(fig12); hold(ax12, 'on'); title('Open loop');
bodeplot(spa_data12_OL_full,"k.");
bodeplot(spa_data12_OL   ,opts1,"r*");%showConfidence(h)
bodeplot(tfest_spa_data12_OL   ,opts1,"b");%showConfidence(h);
bodeplot(spa_OL_from_Tune_and_CL_spa   ,opts1,"g*");
bodeplot(tfest_spa_OL_from_Tune_and_CL_spa   ,opts1,"g-");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","OL from tune and CL"); grid on;

% Closed Loop
fig22 = figure(22);ax22 = axes(fig22); hold(ax22, 'on'); title('Closed loop'); 
bodeplot(spa_data12_CL_full,"k.");
bodeplot(spa_data12_CL   ,opts1,"r*");% showConfidence(h)
bodeplot(tfest_spa_data12_CL   ,opts1,"b");% showConfidence(h);
bodeplot(spa_CL_from_Tune_and_OL_spa   ,opts1,"g*"); 
bodeplot(CL_from_Tune_and_OL_tfest   ,opts1,"g-");
%bodeplot(ssest_data12_CL,opts1,"y-");
bodeplot( CL_from_tune_and_tfest_spa_OL_from_Tune_and_CL_spa,opts1,"y-");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","CL from tune and OL spa","CL from tune and OL tfest"); grid on;

%
% Ymodel = lsim(tfest_spa_data12_OL,sv2_acq,time_acq); 
% E = Ymodel - x_acq_T;
% half=floor(length(time_acq)/2);
% figure(92); hold on; autocorr(E(1:half),NumLags=300);
% %figure; plot(time_acq(1:lags),E)
% % R_XE = xcorr(E,x_acq_T,'coeff'); %  max(E) = 8e+277
% % figure(92); plot( -time_acq(end):Ts:time_acq(end) , R_XE)
Ymodel_tfest_spaOL = lsim(tfest_spa_OL_from_Tune_and_CL_spa,sv2_acq,time_acq); 
E = Ymodel_tfest_spaOL - x_acq_T;
half=floor(length(time_acq)/2);
figure(92); autocorr(E(1:half),NumLags=100);

%% Data 13
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
true_tune_13 = pid(10,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data13_OL = iddata(x_acq_T, sv2_acq, Ts);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data13_CL =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);
%
spa_data13_OL_full = spa(data13_OL, win_size, f_vector_full);
spa_data13_OL = spa(data13_OL, win_size, f_vector_OL);
spa_data13_CL_full = spa(data13_CL, win_size, f_vector_full);
spa_data13_CL = spa(data13_CL, win_size, f_vector_CL);

tfest_spa_data13_OL = tfest(spa_data13_OL,np_OL,'Ts',Ts,tfest_opt_OL)
tfest_spa_data13_CL = tfest(spa_data13_CL,np_CL,'Ts',Ts,tfest_opt_CL)
ssest_data13_CL =  ssest(spa_data13_CL_full,6,ssestOptions)
%
spa_OL_from_Tune_and_CL_spa =  spa_data13_CL_full/(d2d(true_tune_13,Ts)*(1-spa_data13_CL_full));%minreal() 
spa_CL_from_Tune_and_OL_spa = feedback(d2d(true_tune_13,Ts)*spa_data13_OL_full, 1);
CL_from_Tune_and_OL_tfest = feedback(true_tune_13*d2d(tfest_spa_data13_OL,Ts_fpga), 1);

% Open Loop
fig13 = figure(13);ax13 = axes(fig13); hold(ax13, 'on'); title('Open loop');
bodeplot(spa_data13_OL_full,"k.");
bodeplot(spa_data13_OL   ,opts1,"r*");%showConfidence(h)
bodeplot(tfest_spa_data13_OL   ,opts1,"b");%showConfidence(h);
bodeplot(spa_OL_from_Tune_and_CL_spa   ,opts1,"g*");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","OL from tune and CL"); grid on;

% Closed Loop
fig23 = figure(23);ax23 = axes(fig23); hold(ax23, 'on'); title('Closed loop'); 
bodeplot(spa_data13_CL_full,"k.");
bodeplot(spa_data13_CL   ,opts1,"r*");% showConfidence(h)
bodeplot(tfest_spa_data13_CL   ,opts1,"b");% showConfidence(h);
bodeplot(spa_CL_from_Tune_and_OL_spa   ,opts1,"g*"); 
bodeplot(CL_from_Tune_and_OL_tfest   ,opts1,"g-");
bodeplot(ssest_data13_CL,opts1,"y-");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","CL from tune and OL spa","CL from tune and OL tfest"); grid on;

%
Ymodel = lsim(tfest_spa_data13_OL,sv2_acq,time_acq); 
E = Ymodel - x_acq_T;
%figure; plot(time_acq(1:lags),E)
half=floor(length(time_acq)/2);
figure(93); hold on; autocorr(E(1:half),NumLags=300);
% R_XE = xcorr(E,x_acq_T,'coeff'); %  max(E) = 8e+277
% figure(92); plot( -time_acq(end):Ts:time_acq(end) , R_XE)

%% Data 14
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
true_tune_14 = pid(15,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data14_OL = iddata(x_acq_T, sv2_acq, Ts);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data14_CL =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);
%
spa_data14_OL_full = spa(data14_OL, win_size, f_vector_full);
spa_data14_OL = spa(data14_OL, win_size, f_vector_OL);
spa_data14_CL_full = spa(data14_CL, win_size, f_vector_full);
spa_data14_CL = spa(data14_CL, win_size, f_vector_CL);

tfest_spa_data14_OL = tfest(spa_data14_OL,np_OL,'Ts',Ts,tfest_opt_OL)
tfest_spa_data14_CL = tfest(spa_data14_CL,np_CL,'Ts',Ts,tfest_opt_CL)

ssest_data14_CL =  ssest(spa_data14_CL_full,6,ssestOptions)
%
spa_OL_from_Tune_and_CL_spa =  spa_data14_CL_full/(d2d(true_tune_14,Ts)*(1-spa_data14_CL_full));%minreal() 
spa_CL_from_Tune_and_OL_spa = feedback(d2d(true_tune_14,Ts)*spa_data14_OL_full, 1);
CL_from_Tune_and_OL_tfest = feedback(true_tune_14*d2d(tfest_spa_data14_OL,Ts_fpga), 1);

% Open Loop
fig14 = figure(14);ax14 = axes(fig14); hold(ax14, 'on'); title('Open loop');
bodeplot(spa_data14_OL_full,"k.");
bodeplot(spa_data14_OL   ,opts1,"r*");%showConfidence(h)
bodeplot(tfest_spa_data14_OL   ,opts1,"b");%showConfidence(h);
bodeplot(spa_OL_from_Tune_and_CL_spa   ,opts1,"g*");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","OL from tune and CL"); grid on;

% Closed Loop
fig24 = figure(24);ax24 = axes(fig24); hold(ax24, 'on'); title('Closed loop'); 
bodeplot(spa_data14_CL_full,"k.");
bodeplot(spa_data14_CL   ,opts1,"r*");% showConfidence(h)
bodeplot(tfest_spa_data14_CL   ,opts1,"b");% showConfidence(h);
bodeplot(spa_CL_from_Tune_and_OL_spa   ,opts1,"g*"); 
bodeplot(CL_from_Tune_and_OL_tfest   ,opts1,"g-");
bodeplot(ssest_data14_CL,opts1,"y-");
legend("Blackman-Tukey spectral analysis","subset of data to fit model","estimated TF","CL from tune and OL spa","CL from tune and OL tfest"); grid on;

%
Ymodel = lsim(tfest_spa_data14_OL,sv2_acq,time_acq); 
E = Ymodel - x_acq_T;
%figure; plot(time_acq(1:lags),E)
half=floor(length(time_acq)/2);
figure(94); hold on; autocorr(E(1:half),NumLags=300);
% R_XE = xcorr(E,x_acq_T,'coeff'); %  max(E) = 8e+277
% figure(92); plot( -time_acq(end):Ts:time_acq(end) , R_XE)


%% Data Laquila
% % input file
% input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\Targets\';
% file = 'LAquilaReducedScale.tgt'; % load input drv
% LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
% 
% folder_2810 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\28-10-2025\';
% file = 'LAquilla_scl0.18_P5_I3_D0.249.acq'; % load output acq
% %scale = 1;
% LTF_to_TXT_then_load_wSV( file , folder_2810 , 'OutputFolder', folder_2810);
% 
% x_drv_T_0 = x_drv_T_0*1e3; % convert to mm
% x_acq_T = x_acq_T*1e3;
% 
% n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
% data1_CL =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);data13_CL.InputName  = 'x_drv_T_0';data13_CL.OutputName = 'x_acq_T';data13_CL.TimeUnit   = 'seconds';
% 

% %% Closed Loop
% g_data11_CL = spa(data11_CL, win_size);
% g_data12_CL = spa(data12_CL, win_size);
% g_data13_CL = spa(data13_CL, win_size);
% g_data14_CL = spa(data14_CL, win_size);
% 
% % Model training
% nx = 6;
% n4sid_data11_CL = n4sid(data11_CL,nx,'Ts',Ts,n4sidOpt); n4sid_data11_CL.InputName  = data11_CL.InputName;n4sid_data11_CL.OutputName = data11_CL.OutputName;
% n4sid_data12_CL = n4sid(data12_CL,nx,'Ts',Ts,n4sidOpt); n4sid_data12_CL.InputName  = data12_CL.InputName;n4sid_data12_CL.OutputName = data12_CL.OutputName;
% n4sid_data13_CL = n4sid(data13_CL,nx,'Ts',Ts,n4sidOpt);n4sid_data13_CL.InputName  = data13_CL.InputName;n4sid_data13_CL.OutputName = data13_CL.OutputName;
% n4sid_data14_CL = n4sid(data14_CL,nx,'Ts',Ts,n4sidOpt);n4sid_data14_CL.InputName  = data14_CL.InputName;n4sid_data14_CL.OutputName = data14_CL.OutputName;


% %% Figures Open Loop
% fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
% bodeplot(spa_data11_OL   ,opts1, "r*");
% bodeplot(tfest_spa_data11_OL , opts1 , "r-");
% bodeplot(spa_data12_OL   ,opts1, "m*");
% bodeplot(tfest_spa_data12_OL , opts1 , "m-");
% bodeplot(spa_data13_OL   ,opts1, "b*");
% bodeplot(tfest_spa_data13_OL , opts1 , "b-");
% bodeplot(spa_data14_OL   ,opts1, "g*");
% bodeplot(tfest_spa_data14_OL , opts1 , "g-");
% 
% legend(); grid on
% 
% %% Figures Closed Loop
% fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop'); 
% bodeplot(spa_data11_CL   ,opts1, "r*");
% bodeplot(spa30_data11_CL   ,opts1, "ro");
% bodeplot(tfest_spa_data11_CL , opts1 , "r-");
% bodeplot(n4sid_spa_data11_CL , opts1 , "r+-");
% bodeplot(spa_data12_CL   ,opts1, "m*");
% bodeplot(tfest_spa_data12_CL , opts1 , "m-");
% bodeplot(spa_data13_CL   ,opts1, "b*");
% bodeplot(tfest_spa_data13_CL , opts1 , "b-");
% bodeplot(spa_data14_CL   ,opts1, "g*");
% bodeplot(tfest_spa_data14_CL , opts1 , "g-");
% 
% legend(); grid on

% 
% %
% fig11 = figure(11);ax11 = axes(fig11); hold(ax11, 'on'); title('Closed loop');
% h = bodeplot(g_data11_CL   ,opts1, "r*");
% bodeplot(n4sid_data11_CL ,opts1)
% %bodeplot(tfest_data11_CL,opts1)
% %bodeplot(G_PIDF_tracking_data11 , opts1 )
% bodeplot(G_PIDF_true_tune_data11 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
% 
% %
% fig12 = figure(12);ax12 = axes(fig12); hold(ax12, 'on'); title('Closed loop');
% h=bodeplot(g_data12_CL   ,opts1, "m*");
% bodeplot(n4sid_data12_CL ,opts1)
% %bodeplot(G_PIDF_tracking_data12 , opts1 )
% bodeplot(G_PIDF_true_tune_data12 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
% 
% %
% fig13 = figure(13);ax13 = axes(fig13); hold(ax13, 'on'); title('Closed loop');
% h=bodeplot(g_data13_CL   ,opts1, "b*");
% bodeplot(n4sid_data13_CL ,opts1)
% %bodeplot(G_PIDF_tracking_data13 , opts1 )
% bodeplot(G_PIDF_true_tune_data13 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
% 
% %
% fig14 = figure(14);ax14 = axes(fig14); hold(ax14, 'on'); title('Closed loop');
% h=bodeplot(g_data14_CL   ,opts1, "b*");
% bodeplot(n4sid_data14_CL ,opts1)
% %bodeplot(G_PIDF_tracking_data13 , opts1 )
% bodeplot(G_PIDF_true_tune_data14 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;




%% Comparing frequency resolution (window size)
% % spa_data11_OL = spa(data11_OL, win_size, f_vector);
% % spa200_data11_OL = spa(data11_OL,200, f_vector);
% % spa100_data11_OL = spa(data11_OL,100, f_vector);
% % spa30_data11_OL = spa(data11_OL,30, f_vector);
% % spa10_data11_OL = spa(data11_OL,10, f_vector);

% %tfest_spa_data11_OL = tfest(data11_OL,np,'Ts',Ts,tfest_opt);
% % n4sid_spa_data11_OL = n4sid(data11_OL,nx,'Ts',Ts,n4sidOpt);

% % spa_data11_CL = spa(data11_CL, win_size, f_vector);
% % spa200_data11_CL = spa(data11_CL,200, f_vector);
% % spa100_data11_CL = spa(data11_CL,100, f_vector);
% % spa30_data11_CL = spa(data11_CL,30, f_vector);
% % spa10_data11_CL = spa(data11_CL,10, f_vector);

% %tfest_spa_data11_CL = tfest(data11_CL,np,'Ts',Ts,tfest_opt);
% % n4sid_spa_data11_CL = n4sid(data11_CL,nx,'Ts',Ts,n4sidOpt);

% fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
% h = bodeplot(spa_data11_OL   ,opts1,"r.");
% showConfidence(h,3)
% h = bodeplot(spa200_data11_OL,opts1,"c.");
% showConfidence(h,3)
% h = bodeplot(spa100_data11_OL,opts1,"m.");
% showConfidence(h,3)
% h = bodeplot(spa30_data11_OL ,opts1,"b.");
% showConfidence(h,3)
% h = bodeplot(spa10_data11_OL ,opts1,"g.");
% showConfidence(h,3)
% % h = bodeplot(tfest_spa_data11_OL ,opts1);
% % showConfidence(h,3)
% % h = bodeplot(n4sid_spa_data11_OL ,opts1);
% % showConfidence(h,3)
% legend;

% fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop'); 
% h = bodeplot(spa_data11_CL   ,opts1,"r.");
% showConfidence(h)
% h = bodeplot(spa200_data11_CL,opts1,"c.");
% showConfidence(h)
% h = bodeplot(spa100_data11_CL,opts1,"m.");
% showConfidence(h)
% h = bodeplot(spa30_data11_CL ,opts1,"b.");
% showConfidence(h)
% h = bodeplot(spa10_data11_CL ,opts1,"g.");
% showConfidence(h)
% % h = bodeplot(tfest_spa_data11_CL ,opts1);
% % showConfidence(h)
% % h = bodeplot(n4sid_spa_data11_CL ,opts1);
% % showConfidence(h)
% legend;
