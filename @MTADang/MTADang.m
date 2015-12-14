classdef MTADang < MTAData
%MTADang(path,filename,data,sampleRate,model,type,ext)
%
%  MTADang is a subclass of MTAData. 
%
%  Current Data Types: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%
%    third dimension:    marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%    
%    fourth dimension:   spherical coordinates in R3 (theta,phi,r)
%
%
%    Indexing Example:
%       Pitch of 2 markers for all time
%       spine_pitch = ang(:,'spine_middle','spine_upper',2);
%
%       Selected periods for Inter Marker distace
%       spine_pitch = ang([1,300;400,1000],'spine_upper','head_back',3);
%
%  See also MTAData

    properties 
        model
    end

    properties(Transient=true)
        data        % data
    end
    
    methods
        function Data = MTADang(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','ang',[],'angles','a'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
            Data.model = model;
        end
    end
    
end

