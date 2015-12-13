function Data = updatePath(Data,path)
%Data = updatePath(Data,path)        
%
% Inputs:
%   path - string: new path to the directory where the Objects data
%                  should be stored
%
% Outputs:
%   Data - MTAData: Original object passed to this function with an
%                   updated path
%
Data.path = path;
end
