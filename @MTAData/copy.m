function DataCopy = copy(Data,varargin)
% Make a copy of a handle object.
% Instantiate new object of the same class.
%
% TODO varargin should contain a property/value pair for setting 
% and modifying properties during the copy.
%

DataCopy = feval(class(Data),[]);
% Copy all non-hidden properties.
p = properties(Data);
    
    if numel(varargin)>2,
        
    end
    for i = 1:length(p)
        if isa(Data.(p{i}),'MTAData'),
            DataCopy.(p{i}) = Data.(p{i}).copy;
        else
            DataCopy.(p{i}) = Data.(p{i});
        end
    end
end
