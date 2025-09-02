function [Spk]=load_ndm_spikes(Spk, varargin)
% function [Spk]=load_ndm_spikes(Spk, varargin)
%
% A super complex matlab function to load from ndm spike files.
% [ clu, res, fet, spk ]

[MAP_UNT_IND, MAP_GRP_IND, MAP_CLU_IND] = deal(1,2,3);

% >>> DEFARGS >>> -------------------------------------------------------------

defargs = struct(    'units', 'all',                                        ...
                 'incld_fet', false,                                        ...
                 'incld_swf', false,                                        ...
                   'denoise', true                                          ...
);
[units, incld_fet, incld_swf, denoise] =                                    ...
    DefaultArgs(varargin, defargs, '--struct');
% <<< DEFARGS <<< -------------------------------------------------------------

Par = Spk.parent.par;

if ischar(units) && strcmp(units,'all');
    units = Spk.map(:,1)';
elseif ischar(units)
    units = get_unit_set( Spk, Spk.parent, units);
end            

% >>> GET the channel groups contaning the units >>> --------------------------
cluster_groups =                                                            ...
    unique(                                                                 ...
        Spk.map( ismember(Spk.map(:, MAP_UNT_IND), units),                  ...
                MAP_GRP_IND)                                                ...
    );
% <<< GET the channel groups contaning the units <<< --------------------------

% >>> FIND the max count for channels, samples, and features >>> --------------
max_channels = max(                                                         ...
    cellfun(@numel, Par.ElecGp(cluster_groups)));
max_samples  = max(                                                         ...
    [Par.SpkGrps(cluster_groups).nSamples]);
max_features = max(                                                         ...
    arrayfun(@(S) S.nFeatures, Par.SpkGrps(cluster_groups))                 ...
    .* arrayfun(@(S) numel(S.Channels), Par.SpkGrps(cluster_groups)))+1;
% <<< FIND the max count for channels, samples, and features <<< --------------

% >>> LOAD the Spk fields >>> -------------------------------------------------
for cluster_group=cluster_groups(:)'
    load_ndm_clusters(Spk, cluster_group, units, denoise);
    load_ndm_res(Spk, cluster_group);
    if incld_fet && ~isempty(Spk.locs)
        load_ndm_fet( Spk, cluster_group, max_features);
    end
    if incld_swf && ~isempty(Spk.locs)
        load_ndm_swf( Spk, cluster_group, max_channels, max_samples);
    end
    spk.locs = [];
end
% <<< LOAD the Spk fields <<< -------------------------------------------------
   
% >>> SORT the spikes by time >>> ---------------------------------------------
[Spk.res,si] = sort(Spk.res);
Spk.clu = Spk.clu(si);
if incld_fet
    Spk.fet = Spk.fet(si,:);
end
if incld_swf
    Spk.spk = Spk.spk(si,:,:);
end
% <<< SORT the spikes by time <<< ---------------------------------------------

return
