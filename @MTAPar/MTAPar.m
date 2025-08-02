classdef MTAPar < hgsetget
%MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
%
%  MTAPar contains parameters for the recording session.
%
%  
%  Fields:
%
%  Example:
%    
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
%
%    ext:        string, A short, unique string which will be the primary
%                        file identifier
%

    properties 

        %path - string: the directory where the object's data is stored
        path        

        %filename - string: the file name of the .mat file which contains the objects data
        filename

    end
    
    properties( Transient=true )
        %parent - MTASession/MTATrial
        parent = [];
        
        %data - numericArray: (TimeIndex,Feature,Dimension)
        data
        
        treatmentRecord = {};        

    end
    
    methods
        function Par = MTAPar( parent )
            Par.filename = [parent.name, '.xml' ];
            Par.path = parent.spath;
            Par.parent = parent;
        end
    end

end