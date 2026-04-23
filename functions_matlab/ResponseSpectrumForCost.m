function picos_ddx_m = ResponseSpectrumForCost( ref_accel)
    m    = 1;
    zeta = 0.05;

    f_i      = 0.1;
    f_n      = 20;
    n_points = 2e3;
    f_vector = logspace(log10(f_i), log10(f_n), n_points);

    picos_ddx_m = zeros(n_points, 1);

    for i = 1:n_points
        k = m * (2*pi*f_vector(i))^2;
        c = zeta * 2 * m * 2*pi*f_vector(i);

        % Discretize and use filter instead of lsim
        sys_d    = c2d(tf([c k], [m c k]), 0.005, 'zoh');
        [b, a]   = tfdata(sys_d, 'v');
        ddx_m    = filter(b, a, ref_accel);

        picos_ddx_m(i) = max(abs(ddx_m));
    end
end
