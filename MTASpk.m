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
    end
    properties
        %map - matrix: Cluster/electrode mapping
        map
     
        %type - string: Data type
        type
    end
    methods
        
        %s.spk(state,clu,sampleRate,squeeze)
        %s.spk(s.stc{'rear'},
        
        function Spk = MTASpk(varargin)
        % Spk = MTASpk(varargin)
        % The general constuctor for MTASpk 
        % Normal usage is to create an Spk object with all values empty 
        % and then use the create function to populate the properties
        
            [Spk.res,Spk.clu,Spk.map,Spk.sampleRate,Spk.fet,Spk.spk,Spk.type] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],'TimePoints'});            
        end
        function Spk = create(Spk,Session,varargin)
        % Spk = create(Spk,Session,varargin)
        % Function used to load spiking data from {clu,res,fet,spk} file
        % and synchronize with the Session
        %   
        %   varargin:
        %     [sampleRate,states,units,loadField] 
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
        %     loadField: cellarry, NOT IMPLEMENTED
        %
        % See also MTAStateCollection for more information on selecting
        % states
        %
            %% Load and resample Res             
            [sampleRate,states,units,mode,loadField] = DefaultArgs(varargin,{1,[],[],'',{}});
            Spk.sampleRate = sampleRate;
            [Res, Clu, Map] = LoadCluRes(fullfile(Session.spath, Session.name));
            Spk.map = Map;
            Res = Res/Session.sampleRate*Spk.sampleRate;
            if Spk.sampleRate~=1
                Res = ceil(Res);
            end
            Session.sync.resample(1);
            [Res, ind] = SelectPeriods(Res,ceil(Session.sync([1,end])*Spk.sampleRate),'d',1,1);
            Clu = Clu(ind);

            switch mode
                case 'NoBurst'                   
            end
            
            %% Select specific units
            if isempty(units),
                cind = true(numel(Res),1);
            else
                cind = find(ismember(Clu,units));
            end            
            
            Res = Res(cind);
            Clu = Clu(cind);
            
            %% Select specific states
            if ~isempty(states);
                [Res,sind] = SelectPeriods(Res,[Session.stc{states,Spk.sampleRate}],'d',1,0);
                Clu = Clu(sind);
            end
            
            Spk.clu = Clu;
            Spk.res = Res;
        end
        
        function Data = subsref(Data,S)
            ni = numel(S);
            if strcmp(S(1).type,'()')&&ni==1,
                if numel(S.subs)==0,
                    Data = Data.res;
                else
                    Data = Data.res(ismember(Data.clu,S.subs{1}));
                end
                return
            end
            Data = builtin('subsref',Data,S);
        end
        
        function Data = clear(Data)
            Data.res = [];
            Data.clu = [];
            Data.sampleRate = [];
        end
        
        function DataCopy = copy(Data)
        % Make a copy of a handle object.
        % Instantiate new object of the same class.
            DataCopy = feval(class(Data),[]);
            % Copy all non-hidden properties.
            p = properties(Data);
            for i = 1:length(p)
                if isa(Data.(p{i}),'MTAData'),
                    DataCopy.(p{i}) = Data.(p{i}).copy;
                else
                    DataCopy.(p{i}) = Data.(p{i});
                end
            end
        end
                
    end
end