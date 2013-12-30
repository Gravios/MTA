classdef MTADxyz < MTAData
%MTADxyz(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
%
%  MTADxyz is a subclass of MTAData. 
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
%    third dimension:    cartesian coordinates in R3 (x,y,z)
%
%
%    Indexing Example:
%       MTADxyz TimeSeries, xy coordinates of 2 markers for all time
%       xy_head = xyz(:,{'head_back','head_front'},[1,2]);
%
%       MTADxyz TimeSeries, z coordinates of 2 markers for specific periods
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