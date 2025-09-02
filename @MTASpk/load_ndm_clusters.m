function load_ndm_clusters(Spk, cluster_group, units, denoise)
% [Clu, nClusters] = load_nda_cluster(FileName)
%
% A simple matlab function to load a .clu file
%
% .clu file contains rows encoded in ascii
%    - first row, cluster count
%    - second to EOF, cluster id
%

NDM_FILE_EXTENSION = 'clu';

MAP_UNT_IND = 1;
MAP_GRP_IND = 2;
MAP_CLU_IND = 3;

% e.g. file.clu.1
fid = fopen(                             ...
    fullfile(Spk.parent.spath,           ...
             [Spk.parent.name,           ...
              '.',NDM_FILE_EXTENSION,'.',    ...
              num2str(cluster_group)]),  ...
            'r');


if fid==-1,
    error( 'MTA:MTASpk:load_clu:FileIOError' );
end

fseek(fid,0,0);
[nClusters, clu] = deal( fscanf(fid, '%d', 1) ,fscanf(fid, '%d'));
fclose(fid);

if denoise
    clu( clu == 0 | clu == 1 ) = 0;
    nClusters = nClusters - 2;
    if nClusters <= 0
        Spk.locs = [];
        return
    end
end

% >>> Get the clusters and their file positions >>> ---------------------------
% The spike map, Spk.map -> [ unit_id, cluster_group, cluster_id ]

% SELECT the cluster_id(s) from the map
submap = Spk.map( (Spk.map(:, MAP_GRP_IND) == cluster_group), :);

assert(                                                ...
    nClusters == size(submap,1) ,                      ...
    'ncluster does not match the mapped cluster count' ...
);


selected_submap_inds = ismember(submap(:, MAP_UNT_IND), units);

selected_unts  = submap(selected_submap_inds, MAP_UNT_IND);
selected_clus  = submap(selected_submap_inds, MAP_CLU_IND);

% REMAP the clus to the units 
[~, cidx] = ismember(clu, selected_clus);
clu(cidx > 0) = selected_unts( cidx(cidx > 0) );
clu(cidx == 0) = 0;

% $$$ % SELECT the clusters
% $$$ [~, cidx] = ismember(clu, selected_unts);
% $$$ clu(cidx == 0) = 0;

Spk.locs = find( cidx ~= 0 );
Spk.clu = cat(1, Spk.clu, clu(Spk.locs));



% <<< Get the clusters and their file positions <<< ---------------------------


