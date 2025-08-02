function DataCopy = copy(Data)
% Make a copy of a handle object.
% Instantiate new object of the same class.
DataCopy = feval(class(Data),[]);
% Copy all non-hidden properties.
p = properties(Data);
for i = 1:length(p)
    if isa(Data.(p{i}),'MTAData'),
        DataCopy.(p{i}) = Data.(p{i}).copy;
    else
        DataCopy.(p{i}) = Data.(p{i});
    end
end

