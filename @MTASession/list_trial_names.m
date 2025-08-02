function trialNames = list_trial_names(Session)
%trialNames = list_trial_names(Session)
%returns a cell array of trial names associated with the Session

trialNames = {};
files = dir(Session.spath);
re = ['\.trl\.'];
trlFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
for file = 1:length(trlFileList),
    points = regexp(trlFileList{file},'[.]');
    fileparts = [1 points+1; points-1 length(trlFileList{file})]';
    if size(fileparts,1)==5,
        trialNames{end+1} = trlFileList{file}(fileparts(3,1):fileparts(3,2));
    end
end
end
