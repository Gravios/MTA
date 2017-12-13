classdef MTADfbr < MTAData
%MTADfbr(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
%
%  MTADfbr is a subclass of MTAData. 
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
%       chan1and2 = fbr(:,[1,2]);
%
%       Selected periods for the 3rd channel
%       chan3per = fbr([1,300;400,1000],3);
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
        function Data = MTADfbr(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','fbr',[],'fbr','f'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
            Data.model = model;
            
        end        
    end
    
end