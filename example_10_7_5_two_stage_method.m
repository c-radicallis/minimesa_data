% Example 10.7.5 - two stage closed loop identification method
% from Van den Hof lecture notes
clear;clc; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-leftwin_size
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
close all;

% bode plot options
opts2=bodeoptions('cstprefs');opts2.PhaseVisible='on';opts2.MagLowerLimMode = 'manual';opts2.MagLowerLim=-15;opts2.XLim={[0.05 pi]};opts2.PhaseWrapping="on"; %opts2.PhaseWrappingBranch=-360;

z=tf('z');
Ts=1;

%%
data_points=2^15;
t = [0:Ts:data_points-1]';
e_vector = [0.15  , 0.8].*wgn(data_points , 1 , 0.1 );
r1 = 0*wgn(data_points , 1 , 0.1 );
r2 = 1*wgn(data_points , 1 , 0.1 );

win_size = data_points/2^7;
frq_res =1/(win_size*Ts); 

np=2;
nz=2;
G0 = 1/(1-1.6*z^-1+0.89*z^-2)
H0 = (1 - 1.56*z^-1 + 1.045*z^-2 -0.3338*z^-3)/(1 - 2.35*z^-1 + 2.09*z^-2 -0.6675*z^-3)

mar= [];

for i = [1.5 , 0.5 ] %2-stage stops working when S0 is close to instable
    C=i*(z^-1 - 0.8*z^-2)
    S0 = 1/(1+G0*C)
    
    % this is to check the effect of feedback intensity vs excitation intensity
    r = r1 + lsim( C , r2 , t );%reference increases in magnitude when controller gain decreases
    mar = [mar ; mean(abs(r))]
    u_r =  lsim(S0 , r ,t);

    for e = e_vector
        v=lsim(H0 , e, t);
        y=lsim(G0*S0 , r , t)+lsim(S0*H0 , e , t);
        u = r - lsim(C , y ,t);

        % figure;hold on;plot(t,r , t,v , t,y);
        % xlim([0 40]);
        % xlabel('Time');
        % ylabel('Signal')
        % legend('r','v','y'); grid on;

        %% step1
        S_est = polyest(r,u ,  [[0 15 0 0 0] 0],'Ts',Ts);%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%

        % figure; hold on;
        % bodeplot(S0,'r', opts1);
        % bodeplot(S_est,'g', opts1);
        % legend; grid on;

        %% step 2
        u_r_est = lsim(S_est , r , t);

        % figure;hold on;
        % plot(t,r,'*-', 'DisplayName', 'reference')
        % plot(t,v, 'DisplayName', 'v (noise) ')
        % plot(t,u , 'DisplayName', 'u (corrupted)');
        % plot(t,u_r , '-' , 'DisplayName', 'u_r (noise-free)');
        % plot(t,u_r_est ,'g--', 'DisplayName', 'u_r^{estimated}');
        % plot(t,y,'--', 'DisplayName', 'y')
        % xlim([0 12]);
        % xlabel('Time');
        % ylabel('Signal')
        % legend; grid on;

        %% step 3
        spa_100 = spa(iddata(y,u,Ts), win_size)
        %spa_10 = spa(iddata(y,u,Ts), win_size_small)
        G_est = tfest(u_r_est , y , np,nz ,'Ts',Ts,'Feedthrough',true)
        G_direct = tfest( u , y , np,nz,'Ts',Ts,'Feedthrough',true)

        figure;hold on;
        bodeplot(G0,'r', opts2);
        bodeplot( -1/C, 'r--',opts2);ch = get(gca,'Children');set(ch(1),'DisplayName','-1/C');
        bodeplot(spa_100 , 'm.', opts2);
        %bodeplot(spa_10 , 'y.', opts2);
        bodeplot( G_direct , 'b-.' , opts2);
        bodeplot(G_est , 'g--', opts2);
        legend; grid on;

        %% RESID
        %opt_resid = residOptions('MaxLag',20);
        % figure;
        % resid(iddata(y,u,Ts),G0,'r',G_direct,'b-.',G_est,'g--')%,opt_resid);%,OL_direct,OL_indirect
        % legend; grid on;

    end
end