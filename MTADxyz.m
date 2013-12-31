classdef MTADxyz < MTAData
%MTADxyz(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
%
%  MTADxyz is a subclass of MTAData, which stores position data from a vicon
%  recording system.
%
%  NOTE: Future versions will allow other systems
%
%  Current Data Type: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%    
%    third dimension:    cartesian coordinates in R3 (x,y,z)
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
    end
    
    properties(Transient=true)
        data        % data
    end
    
    methods
        function Data = MTADxyz(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','pos'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename= [filename '.' ext '.mat'];
                end
            end            
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext);            Data.model = model;
        end
        function Data = create(Data,varargin)
        %Data = create(Data,varargin)
        %not implemented in this version
        end
        function Data = embed(Data,win,overlap)
        %Data = embed(Data,win,overlap)
        %not implemented in this version
        end

    end
    
end