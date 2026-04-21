function stop = recordAndStop(optimValues, state)
    persistent history
    if strcmp(state, 'init')
        history = [];
    end
    stop = false;
    if strcmp(state, 'iter')
        history = [history; optimValues.iteration, optimValues.fval];
        assignin('base', 'opt_hist', history);
    end
end