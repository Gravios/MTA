function tind = find_trial_index(Trials,filebase)
tind = 0;
if ~isempty(Trials) && ~isempty(filebase)    
    for tind = 1:numel(Trials)
        if isstruct(Trials)
            if strcmp([Trials(tind).sessionName,'.',Trials(tind).mazeName,'.',Trials(tind).trialName], filebase)
                return
            end
        else
            if strcmp(Trials{tind}.filebase, filebase)
                return
            end
        end
    end
end