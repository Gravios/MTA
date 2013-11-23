classdef MTADlfp < MTAData
    properties 
        model = [];
    end
    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADlfp(varargin)
            [path,filename,data,sampleRate,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],'TimeSeries','lfp'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),['.' ext]),
                    filename = [filename '.' ext];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,type,ext);
        end        
        function Data = load(Data,Session,varargin)
            lfpSync = Session.sync.periods(Session.lfp.sampleRate);
            [channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},lfpSync([1,end])});
            Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));
            if isempty(channels)
                channels = Par.(gselect{1})(gselect{2}).Channels(gselect{3});
            end
            Data.data = LoadBinary(Data.fpath,channels,...
                                   Par.nChannels,[],[],[],periods)';
        end
        function Data = create(Data,varargin)
        end
        function Data = updatePath(Data,path)
            Data.path = path;
        end
        function Data = filter(Data,win)
        end
        function Data = resample(Data,newSampleRate,varargin)
            [interp_type] = DefaultArgs(varargin,{'linear'});            
        end
        function Data = embed(Data,win,overlap)
        end        
    end
    
end