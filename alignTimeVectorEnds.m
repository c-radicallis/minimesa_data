function [x_drv_cut, t_acq_aligned] = alignTimeVectorEnds(t_drv, x_drv, t_acq, Ts)

    % Ensure column vectors
    t_drv = t_drv(:);
    t_acq = t_acq(:);
    x_drv = x_drv(:);

    % Lengths
    N_drv = length(t_drv);
    N_acq = length(t_acq);

    % --- Missing time calculation ---
    missing_samples = max(0, N_drv - N_acq);
    missing_time = missing_samples * Ts;

    fprintf('Missing acquisition time: %.6f seconds\n', missing_time);

    % --- Align acquisition time vector ---
    T_drv = t_drv(end);
    T_acq = t_acq(end);

    if T_drv > T_acq
        shift = T_drv - T_acq;
        t_acq_aligned = t_acq + shift;
    else
        t_acq_aligned = t_acq;
    end

    % --- Cut driver signal to match acquisition length ---
    if N_drv > N_acq
        x_drv_cut = x_drv(end - N_acq + 1 : end);
    else
        x_drv_cut = x_drv;
    end

end
