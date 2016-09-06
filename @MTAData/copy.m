function DataCopy = copy(Data,varargin)
% Make a copy of a handle object.
% Instantiate new object of the same class.
%
% TODO varargin should contain a property/value pair for setting 
% and modifying properties during the copy.
%

if iscell(Data),

    DataCopy = cell(size(Data));
    for d = 1:numel(Data),
        DataCopy{d} = feval(class(Data{d}),[]);
        % Copy all non-hidden properties.
        p = properties(Data{d});
        if numel(varargin)>2,
             
        end
            
        for i = 1:length(p)
            if isa(Data{d}.(p{i}),'MTAData'),
                DataCopy{d}.(p{i}) = Data{d}.(p{i}).copy;
            else
                DataCopy{d}.(p{i}) = Data{d}.(p{i});
            end
        end
    end
    
    
else
    
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
