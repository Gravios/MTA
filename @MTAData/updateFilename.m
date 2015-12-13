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
Data.filename = filename;
end
