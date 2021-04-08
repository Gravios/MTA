function [subjects] = load_take_files(Session,records)

subjects = cell([1,numel(records)]);

for r = 1:numel(records),
    if ~isempty(records{r}{2}),
        ds = load(fullfile(Session.spath, Session.maze.name,records{r}{2}));
        % CONTAINS vars: 'csv', 'frames', 'timestamps', 'numFrames', 'numExportedFrames', 'sampleRate', 'subjects'
        subjects{r} = ds.subjects;
    end
end


