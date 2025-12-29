function results = twoStageMethod(...
    Kp , fir_np, np_CL , np_OL,  Ts , opts1,  sv2_acq,  ...
    x_drv_T_0, time_drv_0, time_acq, ...
    x_acq_T)

tfest_opt_CL = tfestOptions('EnforceStability',1);

[x_drv_T_0_cut , time_acq_aligned] = alignTimeVectorEnds(time_drv_0 , x_drv_T_0, time_acq , Ts);

% step1
S_est = polyest(x_drv_T_0_cut,sv2_acq ,  [[0 fir_np 0 0 0] 0],'Ts',Ts)%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
figure; hold on;
bodeplot(S_est,'g', opts1);
legend; grid on;

% step 2
u_r_est = lsim(S_est , x_drv_T_0 , time_drv_0);

figure, hold on;
plot(time_acq_aligned , sv2_acq , 'DisplayName', 'sv2_acq');
plot(time_drv_0,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
% u_r_est_cut = lsim(S_est , x_drv_T_0_cut , time_acq_aligned);
% plot(time_acq_aligned,u_r_est_cut ,'r--', 'DisplayName', 'u_r^{est}');
legend; grid on;

% step 3
G_est = tfest(u_r_est(end - length(time_acq) + 1 : end) , x_acq_T , 4 , 'Ts' , Ts)
G_est_to_CL = feedback(Kp*G_est, 1);
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
legend; grid on;

% Pack results
results.S_est = S_est;
results.u_r_est = u_r_est;
results.G_est = G_est;
results.G_est_to_CL = G_est_to_CL;
results.G_direct = G_direct;
results.G_direct_to_CL = G_direct_to_CL;
results.G_CL = G_CL;
results.G_indirect = G_indirect;


end