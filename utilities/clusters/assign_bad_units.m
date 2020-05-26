
function assign_bad_units(Trial,badUnits,varargin)
%function assign_bad_units(Trial, badUnits, overwrite)
%
% varargin:
%    Trial -        MTATrial: Trial contaning path and metadata
%    badUnits - numericArray: array of bad units
%    overwrite -     Logical: overwrite file with badUnits otherwise append

% DEFARGS ------------------------------------------------------------------------------------------

defargs = struct('overwrite', false);
[overwrite] = DefaultArgs(varargin,defargs,'--struct');

%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

if ~ overwrite
    if exist(fullfile(Trial.spath,[Trial.name,'.bad_units.mat']),'file'),
        ds = load(fullfile(Trial.spath,[Trial.name,'.bad_units.mat']));
    else,
        ds.badUnits = [];
    end
    
    badUnits = unique([ds.badUnits(:);badUnits(:)]);
end

save(fullfile(Trial.spath,[Trial.name,'.bad_units.mat']),'badUnits');

% END MAIN -----------------------------------------------------------------------------------------