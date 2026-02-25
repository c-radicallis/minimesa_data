% Example 10.7.5 - from Van den Hof lecture notes
% two stage closed loop identification method modified for
% reference r2 only, and proportional controller

clear;clc; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-leftwin_size
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
close all;

% bode plot options
opts2=bodeoptions('cstprefs');opts2.PhaseVisible='on';opts2.XLim={[0.2 pi]};opts2.PhaseWrapping="on"; %opts2.PhaseWrappingBranch=-360;
opts2.MagLowerLimMode = 'manual';opts2.MagLowerLim=-15;


Ts=1;z=tf('z',Ts);

%%
data_points=2^15;
t = [0:Ts:data_points-1]';
e_vector = [0.1 , 1 ].*wgn(data_points , 1 , 1 ); %[0.05 , 0.8 ]

r2 = wgn(data_points , 1 , 1 );

np=2;
nz=2;
G0 = 1/(1-1.6*z^-1+0.89*z^-2)
H0 = (1 - 1.56*z^-1 + 1.045*z^-2 -0.3338*z^-3)/(1 - 2.35*z^-1 + 2.09*z^-2 -0.6675*z^-3);

ma_u= [];
ma_r2= [];
ma_v=[];
signal_to_noise= [];

for i =[tf(5)]% [0.05*(z^-1-0.8*z^-2)]% , %tf(0.01),  0.5*(z^-1-0.8*z^-2)] %controller 
    C= i;
    S0 = 1/(1+G0*C)
    
    % this is to check the effect of feedback intensity vs excitation intensity
    r =  lsim( C , r2, t );%reference increases in magnitude when controller gain decreases
    
    ma_r2 = [ma_r2 ; mean(abs(r2))];
    u_r =  lsim(S0 , r ,t);

    for e = e_vector
        v=lsim(H0 , e, t);
        ma_v = [ma_v ; mean(abs(v))];
        y=lsim(G0*S0 , r , t)+lsim(S0*H0 , e , t);
        u = r - lsim(C , y ,t);
        ma_u = [ma_u ; mean(abs(u))];
        signal_to_noise = [signal_to_noise; snr(r,v)];

        % figure;hold on;plot(t,r , t,v , t,y);
        % xlim([0 40]);
        % xlabel('Time');
        % ylabel('Signal')
        % legend('r','v','y'); grid on;

        %% step1
        S_est = polyest(r,u ,  [[0 2^6 0 0 0] 0],'Ts',Ts);%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%
        % step 2
        u_r_est = lsim(S_est , r , t);

        % figure; hold on;
        % bodeplot(S0,'r', opts1);
        % bodeplot(S_est,'g', opts1);
        % legend; grid on;

        % figure;hold on;
        % plot(t,r2,'*-', 'DisplayName', 'r2')
        % plot(t,v,'--', 'DisplayName', 'v (noise) ')
        % plot(t,u , 'DisplayName', 'u (corrupted)');
        % plot(t,u_r , '-' , 'DisplayName', 'u_r (noise-free)');
        % plot(t,u_r_est ,'g--', 'DisplayName', 'u_r^{estimated}');
        % plot(t,y,'--', 'DisplayName', 'y')
        % xlim([500 512]);
        % ylim([-4 4]);
        % xlabel('Time');
        % ylabel('Signal')
        % legend; grid on;

        %% step 3
        win_size_1 = 2^5;
        frq_res_1 =1/(win_size_1*Ts); 
        Nonparametric_1 = spa(iddata(y,u,Ts), win_size_1)

        win_size_2 = data_points/4;
        frq_res_2 =1/(win_size_2*Ts); 
        Nonparametric_2 = spa(iddata(y,u,Ts), win_size_2)

        G_2stage = tfest(u_r_est , y , np,nz ,'Ts',Ts,'Feedthrough',true)
        G_direct = tfest( u , y , np,nz,'Ts',Ts,'Feedthrough',true)

        %%
        figure;hold on;
        bodeplot(G0,'r', opts2);
        bodeplot( -1/C, 'r--',opts2);ch = get(gca,'Children');set(ch(1),'DisplayName','-1/C');
        bodeplot(Nonparametric_1 , 'm.', opts2);ch = get(gca,'Children');set(ch(1),'DisplayName',sprintf('Nonparametric (%.0f lags)', win_size_1));
        bodeplot(Nonparametric_2 , 'c.' , opts2);ch = get(gca,'Children');set(ch(1),'DisplayName',sprintf('Nonparametric (%.0f lags)', win_size_2));

        bodeplot( G_direct , 'b-.' , opts2);
        bodeplot(G_2stage , 'g--', opts2);

        legend; grid on;

        %% RESID
        % %opt_resid = residOptions('MaxLag',20);
        % figure;
        % mscohere(u,y)
        % %resid(iddata(y,u,Ts),G0,'r',G_direct,'b-.',G_2stage,'g--')%,opt_resid);%,OL_direct,OL_indirect
        % legend; grid on;

    end
end