% Test for computation of response spectra

% Response Spectra settings
f_i=0.1; %freq inicial
f_n=20;  %freq final
n_points = 1e3;
f_vector = logspace( log10(f_i) , log10(f_n) , n_points);

tic
[picos_ddx_tgt_T , ~] = ResponseSpectrum( time_vector , x_tgt_T , ddx_tgt_T, f_vector , 0);
elapsed = toc;
fprintf('ResponseSpectrum Elapsed time : %.4f seconds\n', elapsed);

tic
picos_ddx_tgt_T_new  = ResponseSpectrumForCost( time_vector , ddx_tgt_T);
elapsed = toc;
fprintf('ResponseSpectrumForCost Elapsed time : %.4f seconds\n', elapsed);

%%
fig8 = figure(8); grid on;xlabel('Frequency (Hz)');ylabel('Acceleration (m/s^2)');title('Acceleration Response Spectra - Fault Normal');xlim([1 30]);%ylim([0 ceil(max(picos_ddx_acq_T_20hz(1:385,1))) ])
% color1 = 'blue';color2 = 'red' ;color3 = '#EDB120'; color4 = 'black';color5='#77AC30';color6='#00fff7';% Define colors for lines 1/3 and 2/4
RGB =  get(groot,"FactoryAxesColorOrder"); H = compose("#%02X%02X%02X",round(RGB*255)); color1 = 'g';color2 = H(1) ;color3 = H(2) ; color4 = H(3);

grid on; hold on;legend()
plot(f_vector, picos_ddx_tgt_T     ,'.' , 'LineWidth' , 2, 'Color', color1, 'DisplayName', 'Target');% - Normal
plot(f_vector, picos_ddx_tgt_T_new,'.', 'LineWidth' , 2, 'Color', color2);% - Normal

figure;plot(f_vector , picos_ddx_tgt_T -picos_ddx_tgt_T_new)