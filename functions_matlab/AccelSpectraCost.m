%% ── Objective function ────────────────────────────────────────────────────
function J = AccelSpectraCost(log_q, OL_200, plant_aug, integrator,sumblk1 ,n_states , picos_ddx_tgt_ForCost , ddx_tgt , time_vector , disp_limit)
    % Recover weights from log-space
    q   = exp(log_q);           % [Q1 Q2 Q3 Q4 Qi]
    Q   = diag(q);
    R   = 1;
    

    % ── Build LQI controller ──────────────────────────────────────────────
    try
        K_lqi = lqi(OL_200, Q, R);
    catch
        J = 1e12;   % lqi failed → penalise heavily
        return
    end

    K  = K_lqi(1:n_states);
    Ki = K_lqi(end);

    controller = ss([], [], [], -[K, Ki]);
    controller.InputName  = [OL_200.StateName; {'xi'}];
    controller.OutputName = {'i_sv'};

    % ── Close the loop ────────────────────────────────────────────────────
    try
        Optimal_CL = connect(plant_aug, controller, integrator, sumblk1, 'x_ref', 'y_xT');
    catch
        J = 1e12;
        return
    end

    % ── Stability guard ───────────────────────────────────────────────────
    if ~isstable(Optimal_CL)
        J = 1e12;
        return
    end

    % ── Simulate & compute tracking error ────────────────────────────────
    try
        ddx_sim = lsim(Optimal_CL, ddx_tgt, time_vector, 'zoh');
    catch
        J = 1e12;
        return
    end

    % Blow-up guard
    if max(abs(ddx_sim)) > disp_limit * 1.5
        J = 1e12;
        return
    end

    picos_ddx_ForCost = ResponseSpectrumForCost(  ddx_sim );
    % mean-square accel tracking error  (minimise this)
    J = mean((picos_ddx_ForCost - picos_ddx_tgt_ForCost).^2);

    fprintf('  J = %.2e  |  Q = [%s]\n', J, num2str(q, '%.2e  '));
end