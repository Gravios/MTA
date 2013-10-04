function consolidate_batch(Trial,batch_tag,varargin)
[verbose,trialName,clean_up_files,cell_reduced] = DefaultArgs(varargin,{0,'all',0,1});

Trial = MTATrial(Trial,[],trialName);
%Trial = MTATrial(Trial,trialName);


    files = dir(Trial.spath.analysis);
    re = [batch_tag '\.\d\d\d\d\d\d\.'];
    batchFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
    if length(batchFileList)<1, return,end
    object = load([Trial.spath.analysis batchFileList{1}]);
    obj_type = fieldnames(object);
    obj_type = obj_type{1};
    props = properties(object.(obj_type));
    object_count = 1;
    numIterations = object.(obj_type).numIterations;
    for i = 2:length(batchFileList),
        if verbose,fprintf('iteration count: %i\n',i),end
        if batchFileList{i},
            object_count = object_count+1;
            nobject = load([Trial.spath.analysis batchFileList{i}]);
            if ~strcmp(obj_type,fieldnames(nobject)),
                error('Objects are not of the same type')
            end

            for np = 1:length(props),    

                if iscell(nobject.(obj_type).(props{np})),
                    if verbose,fprintf('prop cat: %s\n',props{np}),end
                    object.(obj_type).(props{np}) = cat(ndims(nobject.(obj_type).(props{np})),object.(obj_type).(props{np}),nobject.(obj_type).(props{np}));
                elseif isequal(object.(obj_type).(props{np}),nobject.(obj_type).(props{np})),
                    continue,                    
                elseif ~isempty(object.(obj_type).(props{np})),
                    if verbose,fprintf('prop cat: %s\n',props{np}),end
                    object.(obj_type).(props{np}) = cat(ndims(nobject.(obj_type).(props{np})),object.(obj_type).(props{np}),nobject.(obj_type).(props{np}));
                end
            end
            numIterations = numIterations+nobject.(obj_type).numIterations;
        end
% $$$         if clean_up_files,
% $$$             system(['rm ' Trial.spath.analysis batchFileList{i}]);
% $$$         end
    end

object.(obj_type).numIterations = numIterations;


if verbose,fprintf('Check for cell reduction\n'),tic, end

%% Reduce cell contents to value of a single iteration if all
%% iterations have the same values.

if cell_reduced,
    if verbose,fprintf('Cell reduction selected\n'),end

    for np = 1:length(props),    
        if iscell(object.(obj_type).(props{np}))
            dMap = repmat({1},1,5);
            opSize = ones(1,5);
            opSize(1:ndims(size(object.(obj_type).(props{np})))) =size(object.(obj_type).(props{np}));
            iter_dim = find(opSize==object.(obj_type).numIterations);
            if ~isempty(iter_dim),
                for d = 1:opSize,
                    if d~=iter_dim,
                        dMap{d} = 1:opSize(d);
                    end
                end
                eqkey = false(1,numIterations);
                eqkey(1) = true;
                compCell={};
                compCell = object.(obj_type).(props{np})(dMap{1},dMap{2},dMap{3},dMap{4},dMap{5});
                for i = 2:object.(obj_type).numIterations,
                    dMap{iter_dim} = i;
                    eqkey(i) = isequal(compCell,object.(obj_type).(props{np})(dMap{1},dMap{2},dMap{3},dMap{4},dMap{5}));
                end
                if sum(eqkey)==numIterations,
                    if verbose,fprintf('reducing property: %s\n',props{np}),end
                    object.(obj_type).(props{np}) = compCell;
                end
            end
        end
    end

    if verbose,fprintf('Cell reduction complete\n',props{np}),toc,end
end



object.(obj_type).save(Trial,0)


end



% $$$     Objs = {};
% $$$     prop = properties(batch_object_search);
% $$$     for j = 1:length(Trial.(Object_File_Tag)),
% $$$         matchCnt = 0;
% $$$         numProp = length(prop);
% $$$         for g = 1:length(prop),
% $$$             if isequal(Trial.(Object_File_Tag){j}.(prop{g}),getfield(batch_object_search,prop{g}))&~isempty(getfield(batch_object_search,prop{g}))
% $$$                 matchCnt = matchCnt + 1;
% $$$             else 
% $$$                 numProp = numProp-1;
% $$$             end
% $$$         end
% $$$         if matchCnt == numProp,
% $$$             Objs{end+1} = Trial.(Object_File_Tag){j};
% $$$         end
% $$$     end
