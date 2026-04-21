function stop = recordAndStop(optimValues, state, MSE_PIDF)
    persistent history
    if strcmp(state, 'init')
        history = [];
    end
    stop = false;
    if strcmp(state, 'iter')
        history = [history; optimValues.iteration, optimValues.fval];
        assignin('base', 'opt_hist', history);
        % stop = optimValues.fval < MSE_PIDF;
    end
end