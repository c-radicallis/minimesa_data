% SISO System Identification Example Script
% This script demonstrates how to perform SISO system identification
% using input-output data in MATLAB with the System Identification Toolbox.

clear;clc;close all;
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'

% Load input drv
pink_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
addpath(pink_folder);
% out_dir = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\pink_noise\';
file = 'pink_noise_40Hz_T3mm_0.drv';
% ok = copyAndRenameFile(in_file, out_dir, newName);

func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\'
addpath(func_folder);
LTF_to_TXT_then_load( file , 'InputFolder', pink_folder , 'OutputFolder', pink_folder);

%%
% --- Step 2: Create an iddata object ---
data = iddata(y', u', Ts);

data.InputName  = 'u';
data.OutputName = 'y';
data.TimeUnit   = 'seconds';

% --- Step 3: Choose a model structure ---
% Here we estimate a second-order transfer function
% You can also try ARX, ARMAX, OE, BJ, state-space, etc.
numerator_order   = 2;
denominator_order = 2;

tf_model = tfest(data, numerator_order, denominator_order);

% Alternative ARX estimation:
% na = 2; nb = 2; nk = 1;  % orders and input delay
% arx_model = arx(data, [na nb nk]);

% --- Step 4: Validate the model ---
% Compare measured output with model output
y2 = lsim(tf_model, u, T);

figure;
subplot(2,1,1);
plot(T, y, 'b', T, y2, 'r--');
legend('Measured y','Model Output y2');
title('Model Fit');
xlabel('Time (s)');
ylabel('Output');

subplot(2,1,2);
resid(data, tf_model);
title('Residual Analysis');

% --- Step 5: Display model information ---
disp('Estimated Transfer Function Model:');
display(tf_model);

% --- Optional: Export model for control design ---
sys_id = tf_model;  % Use 'sys_id' in your controller design

% End of script
