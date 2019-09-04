function Session = updatePaths(Session,varargin)
% Session = updatePaths(Session)
% Change Session.path & Session.spath to the current
% MTAPath.mat configuration found in the matlab path.
Session.path = load('MTAPaths.mat');
Session.spath = fullfile(Session.path.project, Session.name);
propList = properties(Session);
    for i = 1:numel(propList),
        if isa(Session.(propList{i}),'MTAData')||isa(Session.(propList{i}),'MTAStateCollection')
            Session.(propList{i}).updatePath(Session.spath);
        end
    end
end
