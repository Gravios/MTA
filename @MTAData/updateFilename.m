function Data = updateFilename(Data,filename)
%Data = updateFilename(Data,filename)        
%
% Inputs:
%   filename - string: new filename for the passed object
%
% Outputs:
%   Data - MTAData: Original object passed to this function with an
%                   updated filename
%
    if ~isempty(filename),
        if isa(filename,'MTASession'),
            filename = [filename.filebase '.' Data.ext '.' Data.label '.' Data.key '.mat'];        
        end
        Data.filename = filename;            

    end            

end
