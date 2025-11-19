clear;clc;close all; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-left
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
Ts_fpga= 1/5000;

%
% spa settings
freq_resolution = 0.1;
win_size = 1/(freq_resolution*Ts); %frequency resolution = 2*pi/(win_size*Ts) [rad/s] = 1/(win_size *Ts)[Hz]
f_vector_OL = logspace( log10(2*freq_resolution*2*pi) , log10(40*2*pi) , 10);
f_vector_CL = logspace( log10(2*freq_resolution*2*pi) , log10(40*2*pi) , 20);

% tfest settings
np_OL = 4; %number of poles for tf est
np_CL = 5;
tfest_opt = tfestOptions('InitialCondition','zero');

%n4sid settings
nx=5;
n4sidOpt = n4sidOptions;
n4sidOpt.N4Weight = 'SSARX'; %allows unbiased estimates when using closed loop data
n4sidOpt.Focus = 'simulation';
n4sidOpt.InitialState = 'zero';

% bode plot options
opts1=bodeoptions('cstprefs');opts1.FreqUnits = 'Hz';opts1.XLim={[freq_resolution 100]};opts1.PhaseWrapping="on";opts1.PhaseWrappingBranch=-360;%opts1.Ylim={[-40 10]};

%% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv
x_drv_T_0 = x_drv_T_0*1e3; % convert t  o mm

%%  Data 11
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
true_tune = pid(5,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data11_openloop = iddata(x_acq_T, sv2_acq, Ts);data11_openloop.InputName  = 'sv2_acq';data11_openloop.OutputName = 'x_acq_T';data11_openloop.TimeUnit   = 'seconds';
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data11_closedloop =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);data11_closedloop.InputName  = 'x_drv_T_0';data11_closedloop.OutputName = 'x_acq_T';data11_closedloop.TimeUnit   = 'seconds';
%
spa_data11_openloop = spa(data11_openloop, win_size, f_vector_OL);
spa_data11_closedloop = spa(data11_closedloop, win_size, f_vector_CL);

%tfest_data11_openloop = tfest(data11_openloop,np_OL,'Ts',Ts,tfest_opt);
% n4sid_data11_openloop = n4sid(data11_openloop,nx,'Ts',Ts,n4sidOpt);
% tfest_data11_closedloop = tfest(data11_closedloop,np_CL,'Ts',Ts,tfest_opt);
% n4sid_data11_closedloop = n4sid(data11_closedloop,nx,'Ts',Ts,n4sidOpt);

tfest_spa_data11_openloop = tfest(spa_data11_openloop,np_OL,'Ts',Ts,tfest_opt);
% n4sid_spa_data11_openloop = n4sid(spa_data11_openloop,nx,'Ts',Ts,n4sidOpt);
tfest_spa_data11_closedloop = tfest(spa_data11_closedloop,np_CL,'Ts',Ts,tfest_opt)
n4sid_spa_data11_closedloop = n4sid(spa_data11_closedloop,nx,'Ts',Ts,n4sidOpt);
%mbj=bj(edat,[4 4 2 2 0]);  % Estimate Box-Jenkins polynomial model using time-domain data


fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
h = bodeplot(spa_data11_openloop   ,opts1,"r.");
showConfidence(h)
% h = bodeplot(tfest_data11_openloop   ,opts1,"g");% showConfidence(h,3);
h = bodeplot(tfest_spa_data11_openloop   ,opts1,"b");
showConfidence(h);
% h = bodeplot(n4sid_data11_openloop ,opts1);% showConfidence(h,3)
% h = bodeplot(n4sid_spa_data11_openloop ,opts1);% showConfidence(h,3)
legend;

fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop'); 
h = bodeplot(spa_data11_closedloop   ,opts1,"r.");% showConfidence(h)
% h = bodeplot(tfest_data11_closedloop   ,opts1,"g");% showConfidence(h);
h = bodeplot(tfest_spa_data11_closedloop   ,opts1,"b");% showConfidence(h);
% h = bodeplot(n4sid_data11_closedloop ,opts1);% showConfidence(h)
h = bodeplot(n4sid_spa_data11_closedloop ,opts1);% showConfidence(h)
legend;

figure;pzmap(spa_data11_closedloop);

%% Comparing frequency resolution (window size)
spa_data11_openloop = spa(data11_openloop, win_size, f_vector);
spa200_data11_openloop = spa(data11_openloop,200, f_vector);
spa100_data11_openloop = spa(data11_openloop,100, f_vector);
spa30_data11_openloop = spa(data11_openloop,30, f_vector);
spa10_data11_openloop = spa(data11_openloop,10, f_vector);

%tfest_spa_data11_openloop = tfest(data11_openloop,np,'Ts',Ts,tfest_opt);
% n4sid_spa_data11_openloop = n4sid(data11_openloop,nx,'Ts',Ts,n4sidOpt);

spa_data11_closedloop = spa(data11_closedloop, win_size, f_vector);
spa200_data11_closedloop = spa(data11_closedloop,200, f_vector);
spa100_data11_closedloop = spa(data11_closedloop,100, f_vector);
spa30_data11_closedloop = spa(data11_closedloop,30, f_vector);
spa10_data11_closedloop = spa(data11_closedloop,10, f_vector);

%tfest_spa_data11_closedloop = tfest(data11_closedloop,np,'Ts',Ts,tfest_opt);
% n4sid_spa_data11_closedloop = n4sid(data11_closedloop,nx,'Ts',Ts,n4sidOpt);

fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
h = bodeplot(spa_data11_openloop   ,opts1,"r.");
showConfidence(h,3)
h = bodeplot(spa200_data11_openloop,opts1,"c.");
showConfidence(h,3)
h = bodeplot(spa100_data11_openloop,opts1,"m.");
showConfidence(h,3)
h = bodeplot(spa30_data11_openloop ,opts1,"b.");
showConfidence(h,3)
h = bodeplot(spa10_data11_openloop ,opts1,"g.");
showConfidence(h,3)
% h = bodeplot(tfest_spa_data11_openloop ,opts1);
% showConfidence(h,3)
% h = bodeplot(n4sid_spa_data11_openloop ,opts1);
% showConfidence(h,3)
legend;

fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop'); 
h = bodeplot(spa_data11_closedloop   ,opts1,"r.");
showConfidence(h)
h = bodeplot(spa200_data11_closedloop,opts1,"c.");
showConfidence(h)
h = bodeplot(spa100_data11_closedloop,opts1,"m.");
showConfidence(h)
h = bodeplot(spa30_data11_closedloop ,opts1,"b.");
showConfidence(h)
h = bodeplot(spa10_data11_closedloop ,opts1,"g.");
showConfidence(h)
% h = bodeplot(tfest_spa_data11_closedloop ,opts1);
% showConfidence(h)
% h = bodeplot(n4sid_spa_data11_closedloop ,opts1);
% showConfidence(h)
legend;



%%
% G_open_from_knownTune_n_CL = minreal( G_closed/(Controller*(1-G_closed))) 
% 
% % PIDF tuning % Resampled to 5000 Hz
% fpga_n4sid_data11_openloop = d2d(n4sid_data11_openloop, Ts_fpga);
% tuner_opts = pidtuneOptions('DesignFocus','reference-tracking');
% cutoff_frequency = 5; % Hz
% PIDF  = pidtune(fpga_n4sid_data11_openloop,'PIDF',cutoff_frequency*2*pi,tuner_opts)
% G_PIDF_5Hz = feedback(PIDF*fpga_n4sid_data11_openloop, 1);

% PIDF  = pidtune(fpga_n4sid_data11_openloop,'PIDF',tuner_opts)
% G_PIDF_tracking_data11 = feedback(PIDF*fpga_n4sid_data11_openloop, 1);

% G_PIDF_true_tune_data11 = feedback(true_tune*fpga_n4sid_data11_openloop, 1);

% checking if the identification is valid (see Duarte Valerio slides)
% Y_est = sim(n4sid_data11_openloop, )
% [c_drv,lags_drv] = xcorr(x_drv_T_0,x_acq_T);

%% Data 12
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P7.acq'; % load output acq
true_tune = pid(7,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data12_openloop = iddata(x_acq_T, sv2_acq, Ts);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data12_closedloop =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);
%
spa_data12_openloop = spa(data12_openloop, win_size, f_vector);
tfest_spa_data12_openloop = tfest(spa_data12_openloop,np,'Ts',Ts,tfest_opt)
spa_data12_closedloop = spa(data12_closedloop, win_size, f_vector);
tfest_spa_data12_closedloop = tfest(spa_data12_closedloop,np,'Ts',Ts,tfest_opt)

%% Data 13
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P10.acq'; % load output acq
true_tune = pid(10,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data13_openloop = iddata(x_acq_T, sv2_acq, Ts);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data13_closedloop =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);
%
spa_data13_openloop = spa(data13_openloop, win_size, f_vector);
tfest_spa_data13_openloop = tfest(spa_data13_openloop,np,'Ts',Ts,tfest_opt)
spa_data13_closedloop = spa(data13_closedloop, win_size, f_vector);
tfest_spa_data13_closedloop = tfest(spa_data13_closedloop,np,'Ts',Ts,tfest_opt)

%% Data 14
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P15.acq'; % load output acq
true_tune = pid(15,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);
x_acq_T = x_acq_T*1e3;
sv2_acq = -sv2_acq; %output is inverted because the wiring is fliped

data14_openloop = iddata(x_acq_T, sv2_acq, Ts);
n1 = numel(x_drv_T_0);n2 = numel(x_acq_T);nmin = min(n1, n2);
data14_closedloop =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);
%
spa_data14_openloop = spa(data14_openloop, win_size, f_vector);
tfest_spa_data14_openloop = tfest(spa_data14_openloop,np,'Ts',Ts,tfest_opt)
spa_data14_closedloop = spa(data14_closedloop, win_size, f_vector);
tfest_spa_data14_closedloop = tfest(spa_data14_closedloop,np,'Ts',Ts,tfest_opt)

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
% data1_closedloop =  iddata(x_acq_T(1:nmin), x_drv_T_0(1:nmin), Ts);data13_closedloop.InputName  = 'x_drv_T_0';data13_closedloop.OutputName = 'x_acq_T';data13_closedloop.TimeUnit   = 'seconds';
% 

% %% Closed Loop
% g_data11_closedloop = spa(data11_closedloop, win_size);
% g_data12_closedloop = spa(data12_closedloop, win_size);
% g_data13_closedloop = spa(data13_closedloop, win_size);
% g_data14_closedloop = spa(data14_closedloop, win_size);
% 
% % Model training
% nx = 6;
% n4sid_data11_closedloop = n4sid(data11_closedloop,nx,'Ts',Ts,n4sidOpt); n4sid_data11_closedloop.InputName  = data11_closedloop.InputName;n4sid_data11_closedloop.OutputName = data11_closedloop.OutputName;
% n4sid_data12_closedloop = n4sid(data12_closedloop,nx,'Ts',Ts,n4sidOpt); n4sid_data12_closedloop.InputName  = data12_closedloop.InputName;n4sid_data12_closedloop.OutputName = data12_closedloop.OutputName;
% n4sid_data13_closedloop = n4sid(data13_closedloop,nx,'Ts',Ts,n4sidOpt);n4sid_data13_closedloop.InputName  = data13_closedloop.InputName;n4sid_data13_closedloop.OutputName = data13_closedloop.OutputName;
% n4sid_data14_closedloop = n4sid(data14_closedloop,nx,'Ts',Ts,n4sidOpt);n4sid_data14_closedloop.InputName  = data14_closedloop.InputName;n4sid_data14_closedloop.OutputName = data14_closedloop.OutputName;

%% Figures Open Loop
fig1 = figure(1);ax1 = axes(fig1); hold(ax1, 'on'); title('Open loop');
bodeplot(spa_data11_openloop   ,opts1, "r*");
bodeplot(tfest_spa_data11_openloop , opts1 , "r-");
bodeplot(spa_data12_openloop   ,opts1, "m*");
bodeplot(tfest_spa_data12_openloop , opts1 , "m-");
bodeplot(spa_data13_openloop   ,opts1, "b*");
bodeplot(tfest_spa_data13_openloop , opts1 , "b-");
bodeplot(spa_data14_openloop   ,opts1, "g*");
bodeplot(tfest_spa_data14_openloop , opts1 , "g-");

legend(); grid on

%% Figures Closed Loop
fig2 = figure(2);ax2 = axes(fig2); hold(ax2, 'on'); title('Closed loop'); 
bodeplot(spa_data11_closedloop   ,opts1, "r*");
bodeplot(spa30_data11_closedloop   ,opts1, "ro");
bodeplot(tfest_spa_data11_closedloop , opts1 , "r-");
bodeplot(n4sid_spa_data11_closedloop , opts1 , "r+-");
bodeplot(spa_data12_closedloop   ,opts1, "m*");
bodeplot(tfest_spa_data12_closedloop , opts1 , "m-");
bodeplot(spa_data13_closedloop   ,opts1, "b*");
bodeplot(tfest_spa_data13_closedloop , opts1 , "b-");
bodeplot(spa_data14_closedloop   ,opts1, "g*");
bodeplot(tfest_spa_data14_closedloop , opts1 , "g-");

legend(); grid on

% 
% %
% fig11 = figure(11);ax11 = axes(fig11); hold(ax11, 'on'); title('Closed loop');
% h = bodeplot(g_data11_closedloop   ,opts1, "r*");
% bodeplot(n4sid_data11_closedloop ,opts1)
% %bodeplot(tfest_data11_closedloop,opts1)
% %bodeplot(G_PIDF_tracking_data11 , opts1 )
% bodeplot(G_PIDF_true_tune_data11 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
% 
% %
% fig12 = figure(12);ax12 = axes(fig12); hold(ax12, 'on'); title('Closed loop');
% h=bodeplot(g_data12_closedloop   ,opts1, "m*");
% bodeplot(n4sid_data12_closedloop ,opts1)
% %bodeplot(G_PIDF_tracking_data12 , opts1 )
% bodeplot(G_PIDF_true_tune_data12 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
% 
% %
% fig13 = figure(13);ax13 = axes(fig13); hold(ax13, 'on'); title('Closed loop');
% h=bodeplot(g_data13_closedloop   ,opts1, "b*");
% bodeplot(n4sid_data13_closedloop ,opts1)
% %bodeplot(G_PIDF_tracking_data13 , opts1 )
% bodeplot(G_PIDF_true_tune_data13 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
% 
% %
% fig14 = figure(14);ax14 = axes(fig14); hold(ax14, 'on'); title('Closed loop');
% h=bodeplot(g_data14_closedloop   ,opts1, "b*");
% bodeplot(n4sid_data14_closedloop ,opts1)
% %bodeplot(G_PIDF_tracking_data13 , opts1 )
% bodeplot(G_PIDF_true_tune_data14 , opts1)
% legend("estimated freq response using spectral analysis","n4sid model","expected CL from identified OL and known controller"); grid on;
