function [clu, nClusters] = load_clu(Spk,cluster_group);
% [Clu, nClusters] = load_clu(FileName)
%
% A simple matlab function to load a .clu file
%
% .clu file contains rows encoded in ascii
%    - first row, cluster count
%    - second to EOF, cluster id
%
EXT = 'clu';

fid = fopen(fullfile(Spk.parent.spath,...
                     [Spk.parent.name,'.',EXT,'.',num2str(cluster_group)]),...
            'r');

if fid==-1,    error( 'MTA:MTASpk:load_clu:FileIOError' ); end
[nClusters, clu] = deal( fscanf(fid, '%d', 1) ,fscanf(fid, '%d'));

fclose(fid);