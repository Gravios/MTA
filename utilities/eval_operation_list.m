function sts = eval_operation_list(sts,stsFuncs)
% function eval_operation_list(sts,stsFuncs)
%
%  A function to perform a series of set operations on the states
%  in sts.
%
%  Note: numel(stsFuncs) should equal numel(sts)-1
%
%  argin:
%    sts: cellArray (MTADepoch) - states
%
%    stsFuncs: cellArray: (char) - operators
%

if isempty(stsFuncs),return,end
while ~isempty(stsFuncs),
    switch stsFuncs{1}
        case '&'
            sts{2} = MTADepoch.intersect(sts(1:2));
            sts(1) = [];
        case '^'
            sts{2} = MTADepoch.intersect(sts(1:2));
            sts(1) = [];
        case '+'
            sts{2} = sts{1}+sts{2};
            %sts{2} = MTADepoch.join(sts(1:2));
            sts(1) = [];
        case '|'
            sts{2} = sts{1}+sts{2};
            %sts{2} = MTADepoch.join(sts(1:2));
            sts(1) = [];
        case '-'
            sts{2} = sts{1}-sts{2};
            sts(1) = [];
    end
    stsFuncs(1) = [];
end
end

