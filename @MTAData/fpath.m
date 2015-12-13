function fpath = fpath(Data)
%fpath = fpath(Data)
% Concatenate the path and filename fields to create the full path
% to the file containing the object's data
%
fpath = fullfile(Data.path,Data.filename);
end
