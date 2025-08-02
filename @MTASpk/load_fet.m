function [fet, nFeatures] = load_fet(Spk, cluster_group)
% [fet, nClusters] = load_fet(Spk, cluster_group)
%
% A simple matlab function to load a .clu file
%
% .clu file contains rows encoded in ascii
%    - first row, cluster count
%    - second to EOF, cluster id
%
EXT = 'fet';

Par = Spk.parent.par;

maxFeatureCount = max(arrayfun(@(S) S.nFeatures, Par.SpkGrps) ...
                 .* arrayfun(@(S) numel(S.Channels), Par.SpkGrps))+1;

fid = fopen(fullfile(Spk.parent.spath,...
                     [Spk.parent.name,'.',EXT,'.',num2str(cluster_group)]),...
            'r');

if fid==-1,  error( 'MTA:MTASpk:load_clu:FileIOError' );  end
[nFeatures, fet] = deal( fscanf(fid, '%d', 1) ,fscanf(fid, '%d'));
fclose(fid);
fet = reshape(fet,nFeatures,[])';


if size(fet,2) < maxFeatureCount
    % >>> PAD Spike Waveform Channels >>> -------------------------------------
    fet = cat(                                                              ...
        2,                                                                  ...
        fet,                                                                ...
        zeros(                                                              ...
            [                                                               ...
                size(fet,1),                                                ...
                maxFeatureCount-size(fet,2)                                 ...
            ]                                                               ...
        )                                                                   ...
    );
    % <<< PAD Spike Waveform Channels <<< -------------------------------------
end


