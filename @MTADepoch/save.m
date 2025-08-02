function out = save(Data,varargin)
% function out = save(Data,varargin)
% 
% saves an MTADepoch object under the location in its path property under the
% filename found in its filename property
% 
% Output:
%     out: boolean - true if save was successful
%
[overwrite] = DefaultArgs(varargin,{0});
out = false;

if ~exist(Data.fpath,'file'),
    save( Data.fpath,'Data','-v7.3');
    out = true;
elseif exist(Data.fpath,'file')&&overwrite,
    warning(['Overwriting: ' Data.fpath]);
    out = true;
    save( Data.fpath,'Data','-v7.3');
else
    warning(['File exists: ' Data.fpath, ' - flag the overwrite option  to save']);
end

end
