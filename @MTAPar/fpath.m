function fpath = fpath(Par)
%fpath = fpath(Par)
% Concatenate the path and filename fields to create the full path
% to the file containing the object's data
%
fpath = resolve_path(fullfile(Par.path,Par.filename));
end
