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
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','lfp',[],'lfp','l'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
            Data.filename = [filename '.' ext];
        end        
    end
    
end