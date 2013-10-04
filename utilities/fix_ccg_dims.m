function fix_ccg_dims(Session,trialName,batch_tag)


Session = MTASession('jg05-20120310');

batch_tag = 'ccg.rdpf_nrhp_sur_par1';

files = dir(Session.spath.analysis);
re = [batch_tag '\.\d\d\d\d\d\d\.'];

batchFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
batch_objects = load([Session.spath.analysis batchFileList{1}]);
seed_type = fieldnames(batch_objects);
props = properties(seed_type);

Bccg = batch_objects.Bccg;

Bccg.posind = permute(batch_objects.Bccg.posind,[2,1]);
Bccg.partition_feature = permute(batch_objects.Bccg.partition_feature,[2,1]);
Bccg.partition_boundaries = permute(batch_objects.Bccg.partition_boundaries,[2,1]);

for i = 2:length(batchFileList),
    if batchFileList{i},
        batch_objects = load([Session.spath.analysis batchFileList{i}]);

        batch_objects.Bccg.posind = permute(batch_objects.Bccg.posind,[2,1]);
        batch_objects.Bccg.partition_feature = permute(batch_objects.Bccg.partition_feature,[2,1]);

        if ~isempty(batch_objects.Bccg.partition_boundaries),
        batch_objects.Bccg.partition_boundaries = permute(batch_objects.Bccg.partition_boundaries,[2,1]);
        Bccg.partition_boundaries = cat(2,Bccg.partition_boundaries,batch_objects.Bccg.partition_boundaries);
        end

        Bccg.posind = cat(2,Bccg.posind,batch_objects.Bccg.posind);
        Bccg.partition_feature = cat(2,Bccg.partition_feature,batch_objects.Bccg.partition_feature);

        Bccg.ccg = cat(5,Bccg.ccg,batch_objects.Bccg.ccg);

    end
end

sbccg = size(Bccg.ccg);
Bccg.numIterations = sbccg(end);

save([Session.spath.analysis Session.filebase '.' batch_tag '.mat'],'Bccg','-v7.3');

