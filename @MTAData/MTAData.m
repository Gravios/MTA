classdef MTAData < hgsetget
%MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
%
%  MTAData is a superclass of most MTA data types. It is a container for 
%  general data types, which alows dynamic refrencing and the use of a set
%  of common functions.
%
%  Current Data Types: 
%    TimeSeries, TimePeriods, TimePoints, SpatialMap
%
%  Current Subclasses:
%    MTADang, MTADepoch, MTADlfp, MTADufr, MTADxyz
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   channel/marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%    
%    Nth dimension:      subspace/channel/marker, :', numeric array of
%                              indicies or string corresponding to one of 
%                              the model labels
%
%    Indexing Example:
%       MTADxyz TimeSeries, xy coordinates of 2 markers for all time
%       xy_head = xyz(:,{'head_back','head_front'},[1,2]);
%
%       MTADxyz TimeSeries, z coordinates of 2 markers for specific periods
%       z_head = xyz([1,300;400,1000],'head_front',3);
%
%       MTADang TimeSeries, pitch of 2 markers for all time
%       spine_pitch = ang(:,'spine_middle','spine_upper',2);
%
%  varargin:
%    
%    path:       string, the directory where the object's data is stored
%
%    filename:   string, The file name of the .mat file which contains the 
%                        objects data
%
%    data:       matrix, Data is your data, ... so put your data here
%
%    sampleRate: double, Sampling rate of the associated data
%        
%    syncPeriods: 
%           MTADepoch,    Time in seconds or the indicies indicating where
%                         the data fits in the Session
%           numericArray, The absolute Recording indicies
%
%    syncOrigin: double, Time or index where the data exits in the overall
%                        session
%
%    type:       string, A short string which denotes the type of data held
%                        by the object
%
%    ext:        string, A short, unique string which will be the primary
%                        file identifier
%
%    sampleRate: double, The sampling rate of the associated data
%


    properties 

        %path - string: the directory where the object's data is stored
        path        

        %filename - string: the file name of the .mat file which contains the objects data
        filename

        parent = [];

        %name - object name (i.e. a subject name jg05)
        name

        %label - string: name associated with the type of data
        label
        
        %key - char: single character used for keyboard shortcuts and indexing
        key
        
        %type - string: a short string which denotes the type of data held by the object
        type        
        
        %ext - string: a short file extension
        ext

        %sampleRate - double: Sampling rate of the associated data
        sampleRate  
        
        %syncPeriods - MTADepoch: Time in seconds indicating where the data fits in the Session
        sync
        
        %syncOrigin - double: time of data origin in seconds with respect to the syncPeriods
        origin = 0;
        
        %hash - string: hash modified by functions acting upon MTAData objects
        hash = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz';
    end
    
    properties( Transient=true )
        treatmentRecord = {};

    end
    
    properties (Abstract)
        %data - numericArray: (TimeIndex,Feature,Dimension)
        data        
        
        %model - MTAModel: holds the names of the features for indexing purposes
        model
        
    end
    
    methods
        function Data = MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key)

            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.' label '.' key '.mat'];                
                end
            end            

            Data.filename = filename;
            Data.path = path;
            Data.data = data;
            Data.type = type;
            Data.ext = ext;
            Data.name = name;
            Data.label = label;
            Data.key = key;
            Data.sampleRate = sampleRate;
            Data.sync = syncPeriods;
            Data.origin = syncOrigin;
            Data.update_hash();
        end
    end
    

    methods (Abstract)        
        Data = create(Data,varargin);      
    end

end