classdef MTADlfp < MTAData
%MTADlfp(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
%
%  MTADlfp is a subclass of MTAData. 
%
%  Current Data Type: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   channel, ':', numeric array of indicies
%                                 representing the loaded channels                              
%
%    Indexing Example:
%       All time for 2 channels
%       chan1and2 = lfp(:,[1,2]);
%
%       Selected periods for the 3rd channel
%       chan3per = lfp([1,300;400,1000],3);
%
%
%  See also MTAData
    properties 
        model = [];
    end
    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADlfp(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','lfp'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),['.' ext]),
                    filename = [filename '.' ext];
                end
            end            
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext);
        end        
        function Data = load(Data,Session,varargin)
            %load(Data,Session,varargin)
            %[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
            [channels,gselect,periods] = DefaultArgs(varargin,...
            {[],{'AnatGrps',1,1},[]});
            if isempty(periods),
                if Session.sync.sampleRate~=1,Session.sync.resample(1);end
                periods = round(Session.sync([1,end]).*Session.lfp.sampleRate);
            end
            Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));
            if isempty(channels)
                channels = Par.(gselect{1})(gselect{2}).Channels(gselect{3});
            end
            Data.data = LoadBinary(Data.fpath,channels,...
                                   Par.nChannels,[],[],[],periods)';
            Data.data(Data.data==0)=1;
            Session.resync(Data);
        end
        function Data = create(Data,Session,varargin)
            %create(Data,Session,varargin)
            %[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
            Data.load(Session,varargin{:});
        end
        function phs = phase(Data,varargin)
            [freq_range,n] = DefaultArgs(varargin,{[6,12],2^11});
            tbp = ButFilter(Data.data(:,:),3,freq_range./(Data.sampleRate/2),'bandpass');
            tbp_hilbert = Shilbert(tbp);
            tbp_phase = phase(tbp_hilbert);
            phs = MTADlfp([],[],tbp_phase,Data.sampleRate,Data.sync.copy,Data.origin);
        end
        function Data = embed(Data,win,overlap)
        %Data = embed(Data,win,overlap)
        %not implemented in this version            
        end        
    end
    
end