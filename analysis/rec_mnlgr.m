

Trial = MTATrial('jg05-20120317');
sampleRate = 30;
fet = fet_lgr(Trial,sampleRate);
states = {'rear','mobile','groom','sit'};

iters = 100;
iter_B = cell([iters,1]);
iter_bp_train = cell([iters,1]);
iter_bp_train_N = cell([iters,1]);
iter_bp_test  = cell([iters,1]);
iter_bp_test_N = cell([iters,1]);
iter_new_smat = cell([iters,1]);

pool = parpool(8,'IdleTimeout',1200);
%matlabpool open 6
rstream = RandStream('mt19937ar','Seed',now);
RandStream.setGlobalStream(rstream);

parfor i = 1:iters,
disp(['iter: ' num2str(i)])
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev1');
fet = fet_lgr(Trial,sampleRate);
aper = Trial.stc{'a'};
Trial.stc{'mobile'};
new_smat = [];
B = [];
bp_train = {};
bp_train_N = [];
bp_test  = {};
bp_test_N = [];
for s = 1:numel(states)

    sper = Trial.stc{states{s}};
    rind = randperm(sper.size(1));
    bp_train{s} = sper.copy;
    bpN = ceil(numel(rind)/2);
    bp_train{s}.data = bp_train{s}(sort(rind(1:bpN)),:);
    bp_train_N(s) = sum(diff(bp_train{s}.data,1,2))./bp_train{s}.sampleRate;
    bp_test{s}  = sper.copy;
    bp_test{s}.data = bp_test{s}(sort(rind(bpN+1:end)),:);
    bp_test_N(s) = sum(diff(bp_test{s}.data,1,2))./bp_test{s}.sampleRate;
    %bp_train{s}.cast('TimeSeries');
    %bp_train{s}.resample(fet);
end


pper = bp_train{1}.copy;
pper.data = sort(pper.data);
fper = pper.copy;
fper = resample(fper.cast('TimeSeries'),fet);
fmat = fper.data;
for s = 2:numel(states),    
    bp_train{i}.data = sort(bp_train{s}.data);
    pper = pper+bp_train{s};
    fper = pper.copy;
    fper = resample(fper.cast('TimeSeries'),fet);
    fmat = fmat+fper.data;
end
pper.fillgaps;
pper.clean;

fmat = abs(fmat-5);
fmat(fmat==5) = 0;
fmat = MTADxyz('data',fmat,'sampleRate',fet.sampleRate);


pper = resample(pper.cast('TimeSeries'),fet);
[B] = mnrfit(fet(pper.data==1&nniz(fmat.data),:),fmat(pper.data==1&nniz(fmat.data)),'model','nominal');
% compute lgr score of original feature
d_state = mnrval(B,fet.data);
for s = 1:size(d_state,2);    
    [~,dind] = max(abs(d_state),[],2);
    new_smat(dind==s,s) = 1;
end
    



iter_B(i) = {B};
iter_bp_train(i) = {bp_train};
iter_bp_train_N(i) = {bp_train_N};
iter_bp_test(i)  = {bp_test};
iter_bp_test_N(i) = {bp_test_N};
iter_new_smat(i) = {new_smat};
end


save(fullfile(Trial.spath,[Trial.filebase '-rec_mnlgr-randClassHeir-1.mat']), ...
              'iter_B',...
              'iter_bp_train',...
              'iter_bp_train_N',...
              'iter_bp_test',...
              'iter_bp_test_N',...
              'iter_new_smat')
delete(fullfile(Trial.spath,[Trial.filebase '-rec_mnlgr-randClassHeir-3.mat']));
ds = load(fullfile(Trial.spath,[Trial.filebase '-rec_mnlgr-randClassHeir-1.mat']));



fileList = dir(Trial.spath);
files = {};
for file = fileList',
    if ~isempty(regexp(file.name,'rec_mnlgr')),
        files = cat(1,files,file.name);
    end
end

ds = [];
for file = files',
    ts = load(fullfile(Trial.spath,file{1}));
    if isempty(ds),
        ds = ts;
        fields = fieldnames(ts);
    else
        for f = fields'
            ts.(f{1})(cellfun(@isempty,ts.(f{1}))) =[];
            ds.(f{1}) = cat(1,ds.(f{1}),ts.(f{1}));
        end
    end
end    
    
iters = numel(ds.iter_B);

iter_cmat = nan([iters,numel(states)+1,numel(states)+1]);
iter_tpr = nan([iters,numel(states)+1]);
iter_ppv = nan([iters,numel(states)+1]);

for k = 1:iters,
    if size(ds.iter_B{k},2)<3,continue,end
bp_test = ds.iter_bp_test{k};
new_smat = ds.iter_new_smat{k};

pper = bp_test{1}.copy;
pper.data = sort(pper.data);
for i = 2:numel(states),    
    bp_test{i}.data = sort(bp_test{i}.data);
    pper = pper+bp_test{i};
end

pper.fillgaps;
pper.clean;


LGRsmat = MTADxyz('data',tnew_smat,'sampleRate',fet.sampleRate);
%LGRsmat.data = LGRsmat.data(:,[3,1,2,4]);
smat = MTADxyz('data',~~stc2mat(Trial.stc,fet,states),'sampleRate',fet.sampleRate);



cmat = zeros([numel(states)+1,numel(states)+1]);
for i = 1:numel(states)+1,
    for j = 1:numel(states)+1,
        if i==numel(states)+1 && j==numel(states)+1,
            cmat(i,j) = sum(double(sum(smat(pper,:),2)==0&sum(LGRsmat(pper,:),2)==0));
        elseif i==numel(states)+1,
            cmat(i,j) = sum(double(sum(smat(pper,:),2)==0&LGRsmat(pper,j)==1));
        elseif j==numel(states)+1,
            cmat(i,j) = sum(double(smat(pper,i)==1&sum(LGRsmat(pper,:),2)==0));            
        else
            cmat(i,j) = sum(double(smat(pper,i)==1&LGRsmat(pper,j)==1));
        end
    end
end



iter_cmat(k,:,:) = cmat;
iter_tpr(k,:) = bsxfun(@rdivide,diag(cmat),sum(cmat,1)')';
iter_ppv(k,:) = bsxfun(@rdivide,diag(cmat),sum(cmat,2))';
end

nanmean(iter_tpr)
figure,hist(iter_tpr(:,1),30)
figure,hist(iter_tpr(:,2),30)
figure,hist(iter_tpr(:,3),30)
figure,hist(iter_tpr(:,4),30)

nanmean(iter_ppv)
figure,hist(iter_ppv(:,1),30)
figure,hist(iter_ppv(:,2),30)
figure,hist(iter_ppv(:,3),30)
figure,hist(iter_ppv(:,4),30)


nanmean(iter_tpr)
edgs = linspace(75,100,100);
figure,s = 1;
subplot(numel(states),1,s),bar(edgs,histc(100*iter_tpr(:,1),edgs),'histc');xlim([75,100]);
Lines(nanmean(iter_tpr(:,s)).*100,[],'r');
title([states{s} ' : ' num2str(round(nanmean(iter_tpr(:,s)).*100,1))]);s=s+1;
subplot(numel(states),1,s),bar(edgs,histc(100*iter_tpr(:,2),edgs),'histc');xlim([75,100]);
Lines(nanmean(iter_tpr(:,s)).*100,[],'r');
title([states{s} ' : ' num2str(round(nanmean(iter_tpr(:,s)).*100,1))]);s=s+1;
subplot(numel(states),1,s),bar(edgs,histc(100*iter_tpr(:,3),edgs),'histc');xlim([75,100]);
Lines(nanmean(iter_tpr(:,s)).*100,[],'r');
title([states{s} ' : ' num2str(round(nanmean(iter_tpr(:,s)).*100,1))]);s=s+1;
subplot(numel(states),1,s),bar(edgs,histc(100*iter_tpr(:,4),edgs),'histc');xlim([75,100]);
Lines(nanmean(iter_tpr(:,s)).*100,[],'r');
title([states{s} ' : ' num2str(round(nanmean(iter_tpr(:,s)).*100,1))]);s=s+1;
xlabel('True Positive Rate');




%% Validation on same animal but different session

Trial = MTATrial('jg05-20120317');
lgr_model = 'MTAC_hand_labeled_rev1_LGR';
states = {'rear','walk','turn','groom','sit'};
fet = fet_lgr(Trial,30);
bhv_lgr(Trial,true,states,fet);
fet = fet_lgr(Trial);
[Stc,d_state] = bhv_lgr(Trial,false,[],fet);


% Other Session
Trial = MTATrial('jg05-20120310');
fet = fet_lgr(Trial);
[Stc,d_state] = bhv_lgr(Trial,false,[],fet,lgr_model);

lgr_model = 'test_rec_mnlgr';
fet = fet_lgr(Trial,30);
bhv_lgr(Trial,true,states,fet,lgr_model);

fet = fet_lgr(Trial);
[Stc,d_state] = bhv_lgr(Trial,false,[],fet,lgr_model);

% Another Session
Trial = MTATrial('Ed05-20140529','all','ont');
fet = fet_lgr(Trial);
[Stc,d_state] = bhv_lgr(Trial,false,[],fet,lgr_model);


lgr_model = 'Ed05-test_rec_mnlgr';
fet = fet_lgr(Trial,30);
bhv_lgr(Trial,true,states,fet,lgr_model);

fet = fet_lgr(Trial);
[Stc,d_state] = bhv_lgr(Trial,false,[],fet,lgr_model);
Stc.updateMode('Ed05-test_rec_mnlgr2');
Stc.save;

% Conf Mat
LGRsmat = MTADxyz('data',~~stc2mat(Stc,fet,states),'sampleRate',fet.sampleRate);
smat    = MTADxyz('data',~~stc2mat(Trial.stc,fet,states),'sampleRate',fet.sampleRate);
cmat = zeros([numel(states)+1,numel(states)+1]);
ind = ':';%Trial.stc{'a'};
for i = 1:numel(states)+1,
    for j = 1:numel(states)+1,
        if i==numel(states)+1 && j==numel(states)+1,
            cmat(i,j) = sum(double(sum(smat(ind,:),2)==0&sum(LGRsmat(ind,:),2)==0));
        elseif i==numel(states)+1,
            cmat(i,j) = sum(double(sum(smat(ind,:),2)==0&LGRsmat(ind,j)==1));
        elseif j==numel(states)+1,
            cmat(i,j) = sum(double(smat(ind,i)==1&sum(LGRsmat(ind,:),2)==0));            
        else
            cmat(i,j) = sum(double(smat(ind,i)==1&LGRsmat(ind,j)==1));
        end
    end
end

cmat
sacc = bsxfun(@rdivide,diag(cmat),sum(cmat,2))
sacc = bsxfun(@rdivide,diag(cmat),sum(cmat,1)')
