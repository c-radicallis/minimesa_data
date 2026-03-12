% experiment for candeias
Ts = 0.005;
G = tf(20*2*pi,[1 20*2*pi],'inputname',"e",'outputname',"y");
lvdt = tf(8*2*pi,[1,8*2*pi]); 
cutoff_frequency=30;
C_tune = pidtune( G , 'PI' , cutoff_frequency*2*pi )
C = pid(1,200)

%%
sys = feedback(G*C , lvdt);

sys_perfect_sensor =  feedback(G*C,1 );

t = 0:1/200:10;

sine_freq= 15 ;% pure sine input
r = sin(2*pi*sine_freq*t);
y = lsim(sys , r , t);
y_perfect_sensor = lsim(sys_perfect_sensor , r , t);

figure; hold on;
plot(t , r);
plot(t,y);
plot(t,y_perfect_sensor);
legend;

figure; bodeplot(sys ,opts1);