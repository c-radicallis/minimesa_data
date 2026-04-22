function [ picos_ddx_m ] = ResponseSpectrumForCost( t_vector , ref_accel)

    m=1;%1kg
    zeta=0.05; %damping ratio (%)s
    % Response Spectra settings
    f_i=0.1; %freq inicial
    f_n=20;  %freq final
    n_points = 1e2;
    f_vector = logspace( log10(f_i) , log10(f_n) , n_points);
    picos_ddx_m = zeros( length(f_vector) ,1 );
    s=tf('s');
    for i=1:length(f_vector)    
      
        % from the desired natural frequency, determining stiffness and damping  
        k = m*(2*pi*f_vector(i))^2; %N/m
        c = zeta*2*m*2*pi*f_vector(i); %N/m/s
          
        % simulacao do modelo 
        ddx_m = lsim( (c*s+k)/(m*s^2+c*s+k) , ref_accel , t_vector ,'zoh'); % Funçao de tranferencia de mola massa amortecedor
        picos_ddx_m(i)=max(abs( ddx_m(:,1) ));  

    end
end
