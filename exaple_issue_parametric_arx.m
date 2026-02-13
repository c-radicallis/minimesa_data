% Example of problem with parametric identification from Van den Hof lecture notes
clear;clc; setappdata(0, 'AutoStagger_LRDown_Last', []);   % ensure first figure starts at top-leftwin_size
set(0, 'DefaultFigureCreateFcn', @autoStagger_LRDown_relSize);
close all;

% bode plot options
opts2=bodeoptions('cstprefs');opts2.PhaseVisible='on';%opts2.MagLowerLimMode = 'manual';%opts2.MagLowerLim=-15;
opts2.XLim={[0.05 pi]};opts2.PhaseWrapping="on"; %opts2.PhaseWrappingBranch=-360;
z=tf('z');
Ts=1;

data_points=2048;
t = [0:Ts:data_points-1]';
e = wgn(data_points , 1 , 0.1 );

r2 = 0.*wgn(data_points , 1 , 0.1 );

win_size = data_points/2;
frq_res =1/(win_size*Ts); 

G0 = 0.1*z^-1/(1-0.8*z^-1)
H0 = (1 - 1.56*z^-1 + 1.045*z^-2 -0.3338*z^-3)/(1 - 2.35*z^-1 + 2.09*z^-2 -0.6675*z^-3);

ma_u= [];
ma_r2= [];
ma_v=[];

C= tf(10);
S0 = 1/(1+G0*C);
% this is to check the effect of feedback intensity vs excitation intensity
r =  lsim( C , r2, t );%reference increases in magnitude when controller gain decreases

ma_r2 = [ma_r2 ; mean(abs(r2))];
u_r =  lsim(S0 , r ,t);


v=lsim(H0 , e, t);
ma_v = [ma_v ; mean(abs(v))]
y=lsim(G0*S0 , r , t)+lsim(S0*H0 , e , t);
u = r - lsim(C , y ,t);
ma_u = [ma_u ; mean(abs(u))];

%% step 3
spa_100 = spa(iddata(y,u,Ts), win_size)
%spa_10 = spa(iddata(y,u,Ts), win_size_small)
G_direct = tfest( u , y , 1,1,'Ts',Ts,'Feedthrough',true)

figure;hold on;
bodeplot(G0,'r', opts2);
bodeplot( -1/C, 'r--',opts2);ch = get(gca,'Children');set(ch(1),'DisplayName','-1/C');
bodeplot(spa_100 , 'm.', opts2);
bodeplot( G_direct , 'b-.' , opts2);
legend; grid on;
