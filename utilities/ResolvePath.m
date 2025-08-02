% function FilePath = ResolvePath(FileName, Fix)
% Resolves the full path to input files in the filebase directory.
% either uses FileBase.fullpath file or, tries to fix it via guess-work of
% the path composition (when on foreign file system)
% on the original filesystem (i.e. bach) with flag Fix it finds original
% path and saves it in FileBase.fullpath file.
% Then it returns FilePath which has RootPath appended to it depending on
% the file system. For bach and lab clients it's nothing, while for cin
% cluster it is /gpfs01/sirota/bach
%
function FilePath = ResolvePath(FileName, varargin)

[ Fix ] = DefaultArgs(varargin,{0});

% decide if we are fixing it for the filename or filebase

if isdir(['../' FileName]) | FileExists([FileName '.xml'])
    FileBase = FileName;
    FileName = [FileName '.xml'];
    Ext = '.xml';
    IsFile = 0;
else
    FileBase = FileName(1:strfind(FileName, '.')-1);
    Ext = FileName(strfind(FileName, '.'):end);
    IsFile = 1;
end

if ~FileExists(FileName) | Fix
    [~, FilePath] = system(['readlink -m -n -s ' FileName]);
else
    if IsFile
        FilePath = FileName;
    else
        FilePath = FileBase;
    end
end

