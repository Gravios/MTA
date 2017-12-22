function randBatch(Session,trialName,batch_seed,batch_tag,consolidate_batch,clean_up_files)

Session = MTASession(Session);
%Session = MTATrial(Session,trialName);


if consolidate_batch,
    batch_object_seed = load([Session.path.root 'batch/' batch_seed]);
    seed_type = fieldnames(batch_object_seed);
    props = properties(seed_type);

    files = dir(Session.path.root 'batch/');
    re = ['\.' batch_tag '\.' batch_object_seed.(seed_type).name '\.\d\d\d\d\d\d\.'];
    batchFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
    batch_objects = load([Session.spath.analysis batchFileList{1}]);
    seed_type = fieldnames(batch_objects);
    props = properties(seed_type);



    for i = 2:length(batchFileList),
        if batchFileList{i},
            nobject = load([Session.spath.analysis batchFileList{i}]);
            if ~strcmp(obj_type,fieldnames(nobject)),
                error('Objects are not of the same type')
            end

            for np = 1:length(props),    
                if isequal(object.(obj_type).(props(np)),nobject.(obj_type).(props(np)),
                    continue,
                elseif length(object.(obj_type).(props(np))~=0
                    object.(obj_type).(props(np)) = cat(length(size(nobject.(obj_type).(props(np)))),object.(obj_type).(props(np)),nobject.(obj_type).(props(np)));
                end
            end

        end
        if clean_up_files,
            system(['rm ' Session.spath.analysis batchFileList{i}]);
        end
    end
end

save([Session.spath.analysis Session.filebase '.' batch_tag '.mat'],'ConArgs','-v7.3');
end



% $$$     Objs = {};
% $$$     prop = properties(batch_object_search);
% $$$     for j = 1:length(Session.(Object_File_Tag)),
% $$$         matchCnt = 0;
% $$$         numProp = length(prop);
% $$$         for g = 1:length(prop),
% $$$             if isequal(Session.(Object_File_Tag){j}.(prop{g}),getfield(batch_object_search,prop{g}))&~isempty(getfield(batch_object_search,prop{g}))
% $$$                 matchCnt = matchCnt + 1;
% $$$             else 
% $$$                 numProp = numProp-1;
% $$$             end
% $$$         end
% $$$         if matchCnt == numProp,
% $$$             Objs{end+1} = Session.(Object_File_Tag){j};
% $$$         end
% $$$     end
