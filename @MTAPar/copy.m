function ParCopy = copy(Par,varargin)
% Make a copy of a handle object.
% Instantiate new object of the same class.
%
% TODO varargin should contain a property/value pair for setting 
% and modifying properties during the copy.
%

    
ParCopy = feval(class(Par),Par.parent);
% Copy all non-hidden properties.
p = properties(Par);
    
if ~isempty(varargin)
    assert(numel(varargin)==1,'MTA:MTAPar:copy:WrongNumberOfInputs');
    if strcmp(varargin{1},'empty'),
        p(~cellfun(@isempty,regexp(p,'^data$'))) = [];
    end
end

for i = 1:length(p)
    if isa(Par.(p{i}),'MTAPar'),
        ParCopy.(p{i}) = Par.(p{i}).copy;
    else
        ParCopy.(p{i}) = Par.(p{i});
    end
end

