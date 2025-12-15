% Example 10.7.5 - two stage closed loop identification method 
clear;clc; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-left
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
%close all;

% bode plot options
opts1=bodeoptions('cstprefs');opts1.MagLowerLimMode = 'manual'; opts1.MagLowerLim=-25;opts1.XLim={[0.05 4]};opts1.PhaseWrapping="on";%opts1.PhaseWrappingBranch=-360;
opts2=bodeoptions('cstprefs');opts2.MagLowerLimMode = 'manual';opts2.MagLowerLim=-20;opts2.XLim={[0.05 4]};opts2.PhaseWrapping="on"; %opts2.PhaseWrappingBranch=-360;

z=tf('z');
t = [0:1:2048-1]';
e = 1*wgn(2048 , 1 , 0.1 );
r1 = 0*wgn(2048 , 1 , 0.1 );
r2 = 1*wgn(2048 , 1 , 0.1 );

C=z^-1 - 0.8*z^-2
r = r1 + lsim( C , r2 , t );

G0 = 1/(1-1.6*z^-1+0.89*z^-2)

H0 = (1 - 1.56*z^-1 + 1.045*z^-2 -0.3338*z^-3)/(1 - 2.35*z^-1 + 2.09*z^-2 -0.6675*z^-3)
v=lsim(H0 , e, t);

S0 = 1/(1+G0*C)

y=lsim(G0*S0 , r , t)+lsim(S0*H0 , e , t);

u = r - lsim(C , y ,t);
u_r =  lsim(S0 , r ,t);

% figure;hold on;
% plot(t,v , t,r , t,y);
% legend;

% step1
S_est = polyest(r,u ,  [[0 15 0 0 0] 0])%tfest(r,u,9)%armax(r,u , [9*[1 1 1] 1],opt)%oe(r ,u , [ 8 8 1 ] )%

figure; hold on;
bodeplot(S0,'r', opts1);
bodeplot(S_est,'g', opts1);
legend; grid on;

%%
% step 2
u_r_est = lsim(S_est , r , t);

% step 3
G_est = tfest(u_r_est , y , 2)
G_direct = tfest( u , y , 2)

figure;hold on;
bodeplot(G0,'r', opts2);
bodeplot(G_est , 'g--', opts2);
bodeplot( G_direct , 'b-.' , opts2);
legend; grid on;

figure, hold on;
plot(t,u , 'DisplayName', 'u');
plot(t,u_r , '-' , 'DisplayName', 'u_r');
plot(t,u_r_est ,'g--', 'DisplayName', 'u_r^{est}');
legend; grid on;
