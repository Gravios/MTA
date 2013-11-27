classdef MTADxyz < MTAData

    properties
        model
    end
    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADxyz(varargin)
            [path,filename,data,sampleRate,model,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],'TimeSeries','pos'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename= [filename '.' ext '.mat'];
                end
            end            
            Data = Data@MTAData(path,filename,data,sampleRate,type,ext);
            Data.model = model;
        end
        
        function Data = create(Data,varargin)
            [Session,overwrite] = DefaultArgs(varargin,{{},{}});
        end
        function Data = filter(Data,win)
            ds = Data.size;
            Data.data = reshape(Filter0(win,Data.data),ds);
        end
        function Data = resample(Data,DataObj)
            if DataObj.isempty, DataObj.load; dlen = DataObj.size(1); end
            uind = round(linspace(round(Data.sampleRate/DataObj.sampleRate),Data.size(1),DataObj.size(1)));
            Data.data = Data.data(uind,:,:);
            Data.sampleRate = DataObj.sampleRate;            
        end
        function Data = embed(Data,win,overlap)
        end

    end
    
end