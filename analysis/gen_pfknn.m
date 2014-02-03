function pfns = Generate_pfknn(Trial,varargin)
% Generate_pfknn
% Automatically iterates trough a list of states, for each creating a 
% k-nearest neighbors spatial rate estimation.
%
%   varargin:
%     [states,selectedUnits] 
%       
%       states: collection of states for which placefields are calculated
%           default:{'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'}
%
%       units: unit ids which correspond to the first dimension of the MTASpk
%              objects map field
%           default:[] - all units
%
[states,units] = DefaultArgs(varargin,{{'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'},[]});

pfns = cell(1,numel(states));
parfor s = 1:numel(states)
    pfns{s} = MTAAknnpfs(Trial,units,states{s},1,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end
