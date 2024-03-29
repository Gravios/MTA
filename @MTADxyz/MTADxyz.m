classdef MTADxyz < MTAData
%MTADxyz(path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext)
%
%  MTADxyz is a subclass of MTAData.
%  MTADxyz objects contain position data.
%  MTADxyz objects are limited to a single subject.
%
%  Current Data Type: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   marker, ':', 
%                                numeric array of indicies
%                                string corresponding to one of the model labels
%                                cellArray of marker label strings
%    
%    third dimension:    cartesian coordinates in R3 (x,y,z)
%                        cartesian coordinates in R2 (x,y)
%
%    Indexing Example:
%       xy coordinates of 2 markers for all time
%       xy_head = xyz(:,{'head_back','head_front'},[1,2]);
%
%       Selected periods for z coordinates of a marker
%       z_head = xyz([1,300;400,1000],'head_front',3);
%
%
%  See also MTAData
    
    properties
        model
        poseIndex = [];
    end
    
    properties(Transient=true)
        data        % data
    end
    
    methods
        function Data = MTADxyz(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','pos',[],'position','x'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
            Data.model = model;
        end


    end
    
end

