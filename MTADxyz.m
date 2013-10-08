classdef MTADxyz < MTAData

    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADxyz(path,data,sampleRate,varargin)
            [type,ext] = DefaultArgs(varargin,{'TimeSeries','pos'});
            if ~strcmp(path(end-3:end),'.mat'),
                path = [path '.' ext '.mat'];
            end
            Data = Data@MTAData(path,data,sampleRate,type,ext);
        end
        
        function Data = create(Data,varargin)
            [Session,overwrite] = DefaultArgs(varargin,{{},{}});
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