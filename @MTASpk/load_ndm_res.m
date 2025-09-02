function load_ndm_res(Spk, cluster_group);
% Res = LoadRes(Spk, cluster_group, locs)
%
% A simple matlab function to load a .res file
NDM_FILE_EXT = 'res';
% BUFFER_SIZE = 64*1024;             % tune this window size

filepath = fullfile(Spk.parent.spath,                                       ...
                    [Spk.parent.name,                                       ...
                     '.', NDM_FILE_EXT, '.',                                ...
                     num2str(cluster_group)]);

fid = fopen(filepath,'r');

% >>> SETUP memory map for spike file >>> -------------------------------------

% >>> BLOCK version >>> -------------------------------------------------------
% $$$ vals = nan(size(locs));
% $$$ 
% $$$ k = 1;
% $$$ while k <= numel(locs)
% $$$     % Start a run
% $$$     runStart = locs(k);
% $$$ 
% $$$     j = k;
% $$$     while j < numel(locs) && (locs(j+1) - runStart) < BUFFER_SIZE
% $$$         j = j + 1;
% $$$     end
% $$$     % Read the contiguous block once
% $$$     blockStart = runStart;
% $$$     blockEnd   = locs(j);
% $$$     blockBytes = blockEnd - blockStart + 8;  % +8 for last double (size of one record)
% $$$     fseek(fid, blockStart, 'bof');
% $$$     raw = fread(fid, blockBytes/8, 'double');  % whole block in one read
% $$$ 
% $$$     % Place values for all indices in this run
% $$$     runLocs   = locs(k:j);
% $$$     relIdx    = (runLocs - blockStart)/8 + 1;  % convert to 1-based indices
% $$$     vals(k:j) = raw(relIdx);
% $$$ 
% $$$     k = j + 1;
% $$$ end
% $$$ fclose(fid);
% <<< BLOCK version <<< -------------------------------------------------------

if fid==-1,    error(['MTA:MTASpk:load_res:FileIOError']); end
%[nTimePoints, Res] = deal( fscanf(fid, '%d', 1) ,fscanf(fid, '%d'));
%nTimePoints = fscanf(fid, '%d', 1);

fseek(fid,0,0);
Res = fscanf(fid, '%d');
fclose(fid);

Spk.res = cat(1, Spk.res, Res(Spk.locs));

% <<< SETUP memory map for spike file <<< -------------------------------------

