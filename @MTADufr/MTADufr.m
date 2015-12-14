classdef MTADufr < MTAData

    properties 
        model = [];
    end

    properties(Transient=true)
        data        % data
    end

    methods
        function Data = MTADufr(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','ufr',[],'ufr','u'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
        end
    end
    
    
end
