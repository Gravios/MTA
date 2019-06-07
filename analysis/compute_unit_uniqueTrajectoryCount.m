function [tcount,scount,tper,states] = compute_unit_uniqueTrajectoryCount(Trial,units,varargin)
%function tcount = compute_unit_uniqueTrajectoryCount(Trial,units)
% 
% Compute the number of independent trajectories passing with2 a 

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sampleRate',                    16,                                            ...
                 'ddz',                           [],                                            ...
                 'pft',                           [],                                            ...
                 'states',                        {{'theta-groom-sit','rear','hloc','hpause',    ...
                                                    'lloc','lpause','groom','sit'}},             ...
                 'distThresh',                    250);
[sampleRate,ddz,pft,states,distThresh] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

if isempty(pft),  pft = pfs_2d_theta(Trial,units);                                          end
if isempty(ddz),  ddz = compute_ddz(Trial,units,pft,'sampleRate',sampleRate);               end

xyz = resample(preproc_xyz(Trial,'trb'),sampleRate);
stcm = stc2mat(Trial.stc.copy(),xyz,states);

spk = Trial.load('spk',sampleRate,'',units,'deburst');    

tbinInds = discretize(1:size(xyz,1),1:sampleRate*3:size(xyz,1));

tcount = nan([numel(units),numel(states),1]);
scount = nan([numel(units),numel(states),1]);
tper   = nan([numel(units),numel(states),2]);
for u = 1:numel(units),
    for s = 1:numel(states),    
        mask = stcm(:,s)==s & abs(ddz(:,u))<distThresh;
        res = spk(units(u));
        res(res>size(mask,1)) = [];
        res = res(mask(res));
        if numel(res)>2,
            tper(u,s,:) = res([1,end])'./sampleRate;
            scount(u,s) = numel(nonzeros(diff(unique(tbinInds(res)))-1));
            tcount(u,s) = numel(nonzeros(diff(unique(tbinInds(mask)))-1));
        end
    end
end
