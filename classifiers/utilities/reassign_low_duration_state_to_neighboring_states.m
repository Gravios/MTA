function [stcMatrix] = reassign_low_duration_state_to_neighboring_states(stcMatrix,state,duration)
for rp = state.data'
    wdur = diff(rp);
    if wdur < duration,
        stcMatrix = reassign_period_to_neighboring_states(rp,stcMatrix);
    end
end
