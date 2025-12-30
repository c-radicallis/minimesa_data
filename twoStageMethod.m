function results = twoStageMethod(...
    Kp , fir_np, np_CL , np_OL,  Ts , opts1,  sv2_acq,  ...
    x_drv_T_0, time_drv_0, time_acq, ...
    x_acq_T)

tfest_opt_CL = tfestOptions('EnforceStability',1);

[x_drv_T_0_cut , time_acq_aligned] = alignTimeVectorEnds(time_drv_0 , x_drv_T_0, time_acq , Ts);

% step1
%S_est = polyest(x_drv_T_0_cut,sv2_acq ,  [[0 fir_np 0 0 0] 0],'Ts',Ts)
% nl = idSaturation('LinearInterval',[min(sv2_acq),max(sv2_acq)]);
% nl = saturation
% model = idnlarx([0 fir_np 0], nl);
% data = iddata(sv2_acq, x_drv_T_0_cut, Ts);
% S_est = nlarx(y, u, Ts, model);

% Estimate HW model
data = iddata( sv2_acq , x_drv_T_0_cut,Ts);
T = getTrend(data,0);
data_detrended = detrend(data , T);
% T_sloped = getTrend(OL_data,1);
% OL_data_detrended_sloped = detrend(OL_data , T_sloped);
% plot(OL_data,OL_data_detrended, OL_data_detrended_sloped)
S_est = polyest(data_detrended,  [[0 fir_np 0 0 0] 0],'Ts',Ts)
sat = idSaturation('LinearInterval',[-16.0920,15.6920]);
sat.Free =[0 0]
S_est_nonLin = nlhw(data_detrended, S_est, [], sat)

figure; hold on;
bodeplot(S_est,'g', opts1);
legend; grid on;

% step 2
u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);
u_r_est_nonLin = sim(S_est_nonLin , x_drv_T_0 );
figure, hold on;
plot(time_acq_aligned , sv2_acq , 'DisplayName', 'sv2_acq');
plot(time_drv_0,u_r_est ,'g', 'DisplayName', 'u_r^{est}');
plot(time_drv_0,u_r_est_nonLin ,'r', 'DisplayName', 'u_r^{est_nonLin}');
% u_r_est_cut = lsim(S_est , x_drv_T_0_cut , time_acq_aligned);
% plot(time_acq_aligned,u_r_est_cut ,'r--', 'DisplayName', 'u_r^{est}');
legend; grid on;

% step 3
G_est = tfest(u_r_est(end - length(time_acq) + 1 : end) , x_acq_T , 4 , 'Ts' , Ts)
G_est_to_CL = feedback(Kp*G_est, 1);

G_est_nonLin = tfest(u_r_est_nonLin(end - length(time_acq) + 1 : end) , x_acq_T , 4 , 'Ts' , Ts)
G_est_nonLin_to_CL = feedback(Kp*G_est_nonLin, 1);

G_direct = tfest( sv2_acq , x_acq_T , np_OL , 'Ts' , Ts)
G_direct_to_CL = feedback(Kp*G_direct, 1);
G_CL = tfest(x_drv_T_0_cut , x_acq_T , np_CL , 'Ts' , Ts,tfest_opt_CL)
G_indirect = G_CL/(Kp*(1-G_CL))

figure;hold on;
bodeplot(G_CL,'y', opts1);
bodeplot(G_est_to_CL , 'g', opts1);
bodeplot( G_direct_to_CL , 'b' , opts1);
bodeplot(G_est , 'g--', opts1);
bodeplot( G_direct , 'b--' , opts1);
bodeplot( G_indirect , 'r--' , opts1);

bodeplot(G_est_nonLin_to_CL , 'g*-', opts1);
bodeplot(G_est_nonLin , 'g*--', opts1);
legend; grid on;

% Pack results
results.S_est = S_est;
% results.u_r_est = u_r_est;
results.G_est = G_est;
results.G_est_to_CL = G_est_to_CL;
results.G_direct = G_direct;
results.G_direct_to_CL = G_direct_to_CL;
results.G_CL = G_CL;
results.G_indirect = G_indirect;


end