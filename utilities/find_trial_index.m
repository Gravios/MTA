function tind = find_trial_index(Trials,filebase)
tind = 0;
if ~isempty(Trials) && ~isempty(filebase)
    for tind = 1:numel(Trials)
        if strcmp(Trials{tind}.filebase, filebase)
            return
        end
    end
end