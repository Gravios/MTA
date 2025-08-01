function [Spk]=load_clu_res(Spk, varargin)
% function [Spk]=load_clu_res(Spk, varargin)
%
% A super complex matlab function to load from the .clu and .res files.
% 
% Usage: [Spk]=LoadCluRes(Spk ,cluster_groups, ClusToLoad)
% Enter the FileBase without the extension .clu or .res
% cluster_groups - list of electrode groups to load (load them all by default)
% 
%     Map is a matrix displaying the correspondance between new cluster numbers (first column) and inital
%     shank # (second column) and cluster # (third column)
% 
% 

Par = Spk.parent.par;

% >>> DEFARGS >>> -------------------------------------------------------------
defargs = struct('cluster_groups',  [1:Par.nElecGps],                       ...
                      'incld_fet', false,                                   ...
                      'incld_spk', false,                                   ...
                   'remove_noise', true                                     ...
);
[cluster_groups, incld_fet, incld_spk, remove_noise] = ...
    DefaultArgs(varargin,{[1:Par.nElecGps],false,false,1});
% <<< DEFARGS <<< -------------------------------------------------------------

[clu, res, fet, spk, map] = deal( [], [], [], [], []);
maxClu=0;
filebase = fullfile(Spk.parent.spath, Spk.parent.name);

for cluster_group=cluster_groups(:)'
    % for each ElGp, load clu and res
    if ~exist([filebase '.clu.' num2str(cluster_group)],'file')
        continue;
    end

    tClu = Spk.load_clu(cluster_group);
    tRes = Spk.load_res(cluster_group);
    if incld_fet
        tFet = Spk.load_fet(cluster_group);
    end
    if incld_spk
        tSpk = Spk.load_spk(cluster_group);
    end

    % >>> REMOVE clusters artifact and noise >>> ------------------------------
    if remove_noise
        indx=(tClu>1);
        if sum(indx)==0
            continue;
        end;
        tClu = tClu(indx);
        tRes = tRes(indx);
        if incld_fet
            tFet = tFet(indx,:);
        end
        if incld_spk
            tSpk = tSpk(indx,:,:);
        end
        
    end %if remove_noise
    % <<< REMOVE clusters artifact and noise <<< ------------------------------

    % CREATE vector of initial tClu and renames cluster id
    [ugini, b, tClu]=unique(tClu); % LIST initial unique values of tClu
    tClu = maxClu+tClu;            % RENUMBER tClu  1+maxClu
    uClu = unique(tClu);           % Unique clus of current electrode group

    % CONCATENATE all the tClu and tRes
    clu = [clu;tClu];
    maxClu=max(clu);
    res = [res;tRes];
    if incld_fet
        fet = [fet;tFet];
    end
    if incld_spk
       spk = cat(1,spk,tSpk);
    end
    
    % CREATE a "map" matrix
    shk=zeros(length(uClu),1) + cluster_group; % electrode #
    map=[map; uClu, shk, ugini];
end
   
% SORT the spikes not to have surprizes later -A
[Spk.res,si] = sort(res);
Spk.clu = clu(si);
Spk.map = map;
if incld_fet
    Spk.fet = fet(si,:);
end
if incld_spk
    Spk.spk = spk(si,:,:);
end
