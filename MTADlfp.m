classdef MTADlfp < MTAData
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
            [channels,gselect,periods] = DefaultArgs(varargin,...
            {[],{'AnatGrps',1,1},round(Session.sync([1,end])*Session.lfp.sampleRate)});
            Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));
            if isempty(channels)
                channels = Par.(gselect{1})(gselect{2}).Channels(gselect{3});
            end
            Data.data = LoadBinary(Data.fpath,channels,...
                                   Par.nChannels,[],[],[],periods)';
            Data.data(Data.data==0)=1;
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
            if DataObj.isempty, DataObj.load; dlen = DataObj.size(1); end
            uind = round(linspace(round(Data.sampleRate/DataObj.sampleRate),Data.size(1),DataObj.size(1)));
            Data.data = Data.data(uind,:);
                        if isa(Data.sync,'MTAData'),
                Data.sync.resample(newSampleRate);
                Data.origin = round(Data.origin/Data.sampleRate*newSampleRate); 
            end     
            Data.sampleRate = DataObj.sampleRate;   
        end
        function Data = embed(Data,win,overlap)
        end        
    end
    
end