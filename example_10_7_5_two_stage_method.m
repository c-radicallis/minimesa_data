% Example 10.7.5 - two stage closed loop identification method 
close all; z=tf('z');

t = [0:1:2048-1]';
e = 0.01*wgn(2048 , 1 , 1 );
r1 = wgn(2048 , 1 , 1 );
r2 = wgn(2048 , 1 , 1 );

C=z^-1 - 0.8*z^-2
r = r1 + lsim( C , r2 , t );

G0 = 1/(1-1.6*z^-1+0.89*z^-2)

H0 = (1 - 1.56*z^-1 + 1.045*z^-2 -0.3338*z^-3)/(1 - 2.35*z^-1 + 2.09*z^-2 -0.6675*z^-3)
v=lsim(H0 , e, t);

S0 = (1+G0*C)^-1

y=lsim(G0*S0 , r , t)+lsim(S0*H0 , e , t);

u = r - lsim(C , y ,t);
u_r =  lsim(S0 , r ,t);

% figure;hold on;
% plot(t,v , t,r , t,y);
% legend;

% step1
S_est = tfest(r,u,7)

figure; hold on;
bodeplot(S0);
bodeplot(S_est);


% step 2
u_r_est = lsim(S_est , r , t);

% step 3
G_est = tfest(u_r_est , y , 3)

figure;hold on;
bodeplot(G0);
bodeplot(G_est);

figure, hold on;
plot(t,u , 'DisplayName', 'u');
plot(t,u_r , '--' , 'DisplayName', 'u_r');
plot(t,u_r_est , 'DisplayName', 'u_r_est');
