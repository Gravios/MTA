function [Res, nTimePoints] = load_res(Spk, cluster_group);
% [Res] = LoadRes(FileName)
%
% A simple matlab function to load a .res file
EXT = 'res';

fid = fopen(fullfile(Spk.parent.spath,...
                     [Spk.parent.name,'.',EXT,'.',num2str(cluster_group)]),...
            'r');

if fid==-1,    error(['MTA:MTASpk:load_res:FileIOError']); end
%[nTimePoints, Res] = deal( fscanf(fid, '%d', 1) ,fscanf(fid, '%d'));
%nTimePoints = fscanf(fid, '%d', 1);
Res = fscanf(fid, '%d');
fclose(fid);
