function load_ndm_fet(Spk, cluster_group, max_features)
% [fet, nClusters] = load_fet(Spk, cluster_group)
%
% A simple matlab function to load a .clu file
%
% .clu file contains rows encoded in ascii
%    - first row, cluster count
%    - second to EOF, cluster id
%

Par = Spk.parent.par;

NDM_FILE_EXT = 'fet';
NDM_FET_DTYPE = char(num2str(Par.nBits, 'int%i'));

% number_of_features 
BLOCK_SIZE = numel(Par.SpkGrps(cluster_group).Channels) * Par.SpkGrps(cluster_group).nFeatures + 1;

filepath = fullfile(Spk.parent.spath,...
                    [Spk.parent.name, ...
                     '.', NDM_FILE_EXT, '.', ...
                     num2str(cluster_group)]);

% CHECK for binary version
if not( exist([filepath,'.bin'], 'file'))
    % >>> GENERATE binary from fet file >>> -----------------------------------
    fid = fopen(filepath, 'r');
    if fid==-1,  error( 'MTA:MTASpk:load_ndm_fet:FileIOError' );  end
    fet = fscanf(fid, '%d');
    fet(1) = [];
    fclose(fid);

    fidw = fopen([filepath,'.bin'],'w');
    fwrite(fidw, fet, NDM_FET_DTYPE);      
    fclose(fidw);
    
    % <<< GENERATE binary from fet file <<< -----------------------------------
end

% >>> LOAD spike features from binary >>> ---------------------------------
MM = memmapfile([filepath,'.bin']);
MM.Format = NDM_FET_DTYPE;
MM.Writable = false;
idx = (Spk.locs-1) * BLOCK_SIZE + [1:(BLOCK_SIZE)];
idx = reshape(idx', [], 1);
fet = reshape( MM.Data(idx), BLOCK_SIZE, [])';
% >>> LOAD spike features from binary >>> ---------------------------------    


% >>> PAD Spike Features >>> --------------------------------------------------
if size(fet,2) < max_features
    fet = cat(                                                              ...
        2,                                                                  ...
        fet,                                                                ...
        zeros(                                                              ...
            [                                                               ...
                size(fet,1),                                                ...
                max_features - size(fet,2)                                  ...
            ]                                                               ...
        )                                                                   ...
    );
end
% <<< PAD Spike Features <<< --------------------------------------------------

Spk.fet = cat(1, Spk.fet, double(fet));
