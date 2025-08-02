function Spk = load(Spk,Session,varargin)
% Spk = load(Spk,Session,varargin)
% Function used to load spiking data from {clu,res,fet,spk} file
% and synchronize with the Session
%   
%   varargin:
%     [sampleRate,states,units,mode,loadField] 
%
%     sampleRate: double, the sampleRate in Hertz
%
%     states:     string, state label or expression to select for a
%                         subset of spikes
%
%     units:       array, Cluster identities to select for a subset
%                         of spikes. The default [] will return all
%                         Clusters.
%
%     mode:       string, {'','deburst'}, only keep 1st spk of each
%                         burst.
%
%     loadField: cellarry, NOT IMPLEMENTED
%
% See also MTAStateCollection for more information on selecting
% states
%
Spk.create(Session,varargin{:});

