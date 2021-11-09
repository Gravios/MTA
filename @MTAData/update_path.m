function Data = update_path(Data,path)
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
    if isa(path,'MTASession')
        Data.path = path.spath;
    else
        Data.path = path;
    end
end
