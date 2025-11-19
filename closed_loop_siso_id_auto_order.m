clear;clc;close all; 
addpath 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model'
func_folder  =  'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\';
addpath(func_folder);
Ts = 0.005;
Ts_fpga= 1/5000;

% input file - pink noise 40hz
input_file_folder ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\31-7-2025\tgt and noise drv\';
file = 'pink_noise_40Hz_T3mm_0.drv'; % load input drv
LTF_to_TXT_then_load( file , 'InputFolder', input_file_folder , 'OutputFolder', input_file_folder); % load input drv

r = x_drv_T_0*1e3; % convert t  o mm

%  Data 11
folder_0711 ='C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\minimesa_data\7-11-2025\';
file = 'pink_noise_40Hz_T3mm_0_P5.acq'; % load output acq
true_tune = pid(5,0,0,0.0019455 , Ts_fpga  );
LTF_to_TXT_then_load_wSV( file , folder_0711 , 'OutputFolder', folder_0711);

y = x_acq_T*1e3;
uc = -sv2_acq; %output is inverted because the wiring is fliped


% closed_loop_SISO_ID_auto_order.m
% --------------------------------------------------
% Automatic model-order selection for closed-loop SISO
% User provides:
%   r  - reference signal (column vector)
%   uc - controller output (plant input) (column vector)
%   y  - plant output (column vector)
%   Ts - sample time (scalar)
%
% Notes:
% - The script treats `uc` as the plant input and `y` as the plant output.
% - The reference `r` is used for diagnostics / optional instruments but
%   is NOT used as a plant input (so the identified model corresponds to
%   the true plant from uc -> y).
% - Workflow: split data -> IV ARX baseline -> n4sid sweep (1..Nmax) ->
%   pick best order by validation loss -> pem refinement -> validation.
% - Requires System Identification Toolbox.
% --------------------------------------------------

%% ---------- USER: supply your data here ----------
% Example placeholders (replace with your real data):
% r  = your_reference_vector;   % column vector
% uc = your_controller_output;  % column vector (plant input)
% y  = your_plant_output;       % column vector
% Ts = your_sample_time;        % sample time in seconds

% Check that the user replaced the placeholders
if ~exist('r','var') || ~exist('uc','var') || ~exist('y','var') || ~exist('Ts','var')
    error('Please set r, uc, y, and Ts in the workspace before running this script.');
end

% Ensure column vectors
r = r(:); uc = uc(:); y = y(:);
N = min([length(r), length(uc), length(y)]);
if any([length(r), length(uc), length(y)] ~= N)
    warning('Input vectors have different lengths. Truncating to the shortest (%d samples).', N);
    r  = r(1:N);
    uc = uc(1:N);
    y  = y(1:N);
end

%% Build iddata object (SISO: input = uc, output = y)
data = iddata(y, uc, Ts);
%% Optional: attach the measured reference as an extra signal for convenience
% You can keep the reference in the workspace for plotting or to form
% instruments if you want to run IV methods manually.

%% Split into estimation / validation
valFraction = 0.3;                 % change if desired
Nval = max(1,round(valFraction*N));
Nest = N - Nval;
data_est = data(1:Nest);
data_val = data(Nest+1:end);

%% 1) Baseline: Instrumental-variable ARX (fast, consistent under closed loop)
na = 6; nb = 6; nk = 1;   % quick baseline orders (adjust if you want)
try
    sys_iv = iv4(data_est, [na nb nk]);
    disp('IV (iv4) baseline model created.');
catch ME
    warning('iv4 failed or not available: %s\nSkipping IV baseline.', ME.message);
    sys_iv = [];
end

%% 2) Subspace sweep: try orders 1..Nmax and select by validation loss
maxOrder = 10;                % upper bound for model order search
opt = n4sidOptions;           % default options
% prefer simulation accuracy which often helps in closed-loop
try
    opt.Focus = 'simulation';
    opt.EnforceStability = true;
catch
    % older toolbox versions may not support option names in this exact way
end

lossVals = nan(maxOrder,1);
models = cell(maxOrder,1);

for n = 1:maxOrder
    try
        % Estimate state-space model of order n using n4sid
        ssTemp = n4sid(data_est, n, opt);
        models{n} = ssTemp;
        % Compute validation loss (mean squared prediction error / simulation loss)
        % The 'loss' function returns a scalar measure (smaller is better).
        lossVals(n) = loss(ssTemp, data_val);
        fprintf('Order %2d: loss = %g\n', n, lossVals(n));
    catch ME
        warning('n4sid failed for order %d: %s', n, ME.message);
        lossVals(n) = Inf;
        models{n} = [];
    end
end

% Pick best order (lowest loss) and display a small table
[bestLoss, bestIdx] = min(lossVals);
bestOrder = bestIdx;
fprintf('\nSelected model order = %d (validation loss = %g)\n', bestOrder, bestLoss);

sys_n4_best = models{bestOrder};

%% Optional: show Hankel singular values for the chosen model
try
    hsv = hsvd(sys_n4_best);
    figure; semilogy(1:length(hsv), hsv, 'o-');
    title(sprintf('Hankel singular values (order %d)', bestOrder));
    xlabel('state index'); ylabel('HSV'); grid on;
catch
    % ignore if hsvd not available
end

%% 3) PEM refinement starting from the selected subspace model
try
    sys_pem = pem(data_est, sys_n4_best);
    fprintf('PEM refinement finished.\n');
catch ME
    warning('PEM refinement failed: %s\nUsing the subspace model as final model.', ME.message);
    sys_pem = sys_n4_best;
end

%% 4) Compare models on validation data
figure('Name','Model comparison (validation)');
% Plot the three models if available: IV, N4SID best, PEM
modelsToCompare = {};
modelsNames = {};
if ~isempty(sys_iv)
    modelsToCompare{end+1} = sys_iv; modelsNames{end+1} = 'IV (iv4)';
end
if ~isempty(sys_n4_best)
    modelsToCompare{end+1} = sys_n4_best; modelsNames{end+1} = sprintf('N4SID (order %d)',bestOrder);
end
if ~isempty(sys_pem)
    modelsToCompare{end+1} = sys_pem; modelsNames{end+1} = 'PEM (refined)';
end

if isempty(modelsToCompare)
    error('No models available to compare.');
else
    compare(data_val, modelsToCompare{:});
    legend(modelsNames,'Location','Best');
end

%% 5) Residual analysis on the final chosen model (PEM if exists)
figure('Name','Residuals (final model)');
resid(data_val, sys_pem);

%% 6) Save final models and summary
save('identified_models_closedloop.mat', 'sys_iv', 'sys_n4_best', 'sys_pem', 'bestOrder', 'lossVals');

fprintf('\nSaved models to identified_models_closedloop.mat\n');

%% 7) Further suggestions / next steps (printed)
fprintf([
    '\nNext steps and tips:\\n', ...
    '- Inspect residuals for whiteness and lack of correlation to measured inputs.\\n', ...
    "- If residuals show correlation, consider richer noise models (Box-Jenkins) or include measured disturbances.\\n", ...
    "- If you have a model of the controller (L), you can transform to an open-loop identification problem for improved results.\\n", ...
    "- To use the reference r as an instrument in IV methods, you can form instruments from delayed r and uc signals. Let me know if you want that added.\\n" ...
]);

% End of script
