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
        end
        function Data = resample(Data,newSampleRate,varargin)
            [interp_type] = DefaultArgs(varargin,{'linear'});            
        end
        function Data = embed(Data,win,overlap)
        end

    end
    
end