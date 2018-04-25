classdef MTASpk < hgsetget
% MTASpk(varargin)
% Data structure to organize timing and features of neuron action potentials.
%
%   Note: The properties (clu,res,fet,spk) should all have the same length
%         along the first dimension.
%   Note: MTASpk does not have a sync property. At the moment, anytime
%         spiking information is needed, it should be generated from the
%         orginal files and synced to the current MTATrial or MTASession.
%   Note: MTASpk.create only loads clu res and map at the moment.
%   Note: MTASpk.create is only compatible with the Neuralynx system.
%
%   varargin:
%     [res,clu,map,sampleRate,fet,spk,type] 
%
%     res:    columnArray(spkIndex,Timestamp), spike timing index in the 
%                                              specified sampleRate 
%
%     clu:    columnArray(spkIndex,ClusterId), spike cluster identity for 
%                                              each spike
%
%     map:    matrix     (ClusterId,ElectrodeGroup,ElectrodeClusterId)
%                                              Cluster/electrode mapping
%
%     sampleRate: double, the sampling rate to which the spike indicies
%                         correspond 
%
%     fet:    matrix     (spkIndex,Feature),   First couple of spike waveform
%                                              principle components for each
%                                              electrode
%
%     spk:    matrix     (spkIndex,spkWaveform), Filtered spike waveforms
%
%     type:   string, Data type

    properties(Transient=true)
        %clu - columnArray: spike cluster identity for each spike
        clu
        
        %res - columnArray: spike timing index in the specified sampleRate 
        res
        
        %fet - matrix: spike waveform features
        fet

        %spk - matrix: Filtered spike waveforms
        spk
        
        %sampleRate - double: the sampling rate to which the spike indicies correspond
        sampleRate    
        
        %hash - string: hash modified by functions acting upon MTAData objects        
        hash = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz';
    end
    
    properties
        %map - matrix: Cluster/electrode mapping
        map
        
        %per - matrix: (nx2) numeric: computational periods of major shift events
        per
        
        %perInd - matrix: (uxn) logical: included periods
        perInd
        
        %type - string: Data type
        type
    end
    
    methods
        
        function Spk = MTASpk(varargin)
        % Spk = MTASpk(varargin)
        % The general constuctor for MTASpk 
        % Normal usage is to create an Spk object with all values empty 
        % and then use the create function to populate the properties
            [Spk.res,Spk.clu,Spk.map,Spk.sampleRate,Spk.fet,Spk.spk,Spk.type] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],'TimePoints'});            
        end
                
    end%methods
end%classdef