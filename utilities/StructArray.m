
% Out = StructArray(In,dim)
%
% Flips over a struct array.  So if you have a structure where each field
% is an array, it will turn it into a structure array where each member of
% the array is a structure.  And vice versa
%
% where struct -> structArray, dim is the dimension of the array
% you wish to correspond to in the structArray

function Out = StructArray(In,varargin)
[dim] = DefaultArgs(varargin,{[]});

fields = fieldnames(In(1));
Out = InitStructArray(fields,0);
if size(In)==[1 1]
    % we are dealing with a single structure whose fields are arrays

    if isempty(dim),
        Out = Out(ones(length(getfield(In,fields{1})),1));    
    else
        sa = size(getfield(In,fields{1}));
        Out = Out(ones(sa(dim),1));    
    end
    for i=1:length(fields)
        f = fields{i};
        a = getfield(In, f);
        sa = size(a);
        if isempty(dim),
            [dimSize,dimInd] = max(sa);
        else
            dimSize = sa(dim);
            dimInd = dim;
        end
        field_ind = {};
        for k = 1:length(sa),
            field_ind{k} = 1:sa(k);
        end
        field_ind{dimInd} = 1;

        for j=1:dimSize,
            field_ind{dimInd} = j;
            Out = setfield(Out, {j}, f, sq(getfield(In, f, field_ind)));
        end
    end
else
    % this is easier - turn a struct array into a single struct

    for i=1:length(fields)
        f = fields{i};
        for j=1:size(In,1)
            for k=1:size(In,2)
                Val = getfield(In, {j,k}, f);
                if ~isempty(Val)
                    Out = setfield(Out, f, {j,k}, Val);
                end
            end
        end
    end
end

