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
    
if ~isempty(varargin)
    assert(numel(varargin)==1,'MTA:MTAData:copy:WrongNumberOfInputs');
    if strcmp(varargin{1},'empty'),
        p(~cellfun(@isempty,regexp(p,'^data$'))) = [];
    end
end


for i = 1:length(p)
    if isa(Data.(p{i}),'MTAData'),
        DataCopy.(p{i}) = Data.(p{i}).copy;
    else
        DataCopy.(p{i}) = Data.(p{i});
    end
end

