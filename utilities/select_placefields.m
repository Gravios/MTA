function units = select_placefields(Trial,varargin)
%function units = select_placefields(Trial,varargin)
%
% select quality placefield units 
%
%  varargin:
%    minSpkCnt - Numeric: default(30)
%    overwrite - Logical: default(false)
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('minSpkCnt',                    30,                                             ...
                 'overwrite',                    false                                           ...
);
[minSpkCnt,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% TAG creation -------------------------------------------------------------------------------------
% ID Vars - create hash tag

% $$$ tag = 'none'
% $$$ if isempty(tag),
% $$$     tag = DataHash(struct());
% $$$ end

filename = fullfile(Trial.spath,[Trial.filebase,'.select_placefields.mat']);

%---------------------------------------------------------------------------------------------------





% MAIN ---------------------------------------------------------------------------------------------

if exist(filename,'file') && ~overwrite
    load(filename,'units');
else,
    spk = Trial.spk.copy();
    Trial.load('nq');    
    pft = pfs_2d_theta(Trial,[],false,true,1);
    mrt = pft.maxRate(spk.map(:,1));
    units = select_units(Trial,'pyr');
    units = units(mrt(units)>0.5);
    units = units(Trial.nq.Refrac(units)<0.001);
    spk.create(Trial,[],'theta-groom-sit',units,'deburst');
    spkCnt = accumarray(spk.clu,ones([numel(spk.clu),1]),[size(spk.map,1),1],@sum);
    units = units(spkCnt(units)>minSpkCnt);
    save(filename,'units');
end
%---------------------------------------------------------------------------------------------------