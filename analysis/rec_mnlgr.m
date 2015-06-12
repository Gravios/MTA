

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

if exist('s','var'),delete(s);end
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
parfor i = 1:iters,
disp(['iter: ' num2str(i)])
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev1');
fet = fet_lgr(Trial,sampleRate);
aper = Trial.stc{'a'};
Trial.stc{'mobile'};
new_smat = [];
B = {};
bp_train = {};
bp_train_N = [];
bp_test  = {};
bp_test_N = [];
for s = 1:numel(states)
    disp(['fiting state: ' states{s}])
    atemper = aper.copy;
    sper = Trial.stc{states{s}};
    rind = randperm(sper.size(1));
    bp_train{s} = sper.copy;
    bp_train_N(s) = ceil(numel(rind)/2);
    bp_train{s}.data = bp_train{s}(sort(rind(1:bp_train_N(s))),:);
    bp_test{s}  = sper.copy;
    bp_test_N(s) = numel(rind)-bp_train_N(s);
    bp_test{s}.data = bp_test{s}(sort(rind(bp_train_N(s)+1:end)),:);
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

fmat = MTADxyz('data',fmat,'sampleRate',fet.sampleRate);


pper = resample(pper.cast('TimeSeries'),fet);
[B{s}] = mnrfit(fet(pper.data==1&nniz(fmat.data),:),fmat(pper.data==1&nniz(fmat.data)),'model','nominal');
    % compute lgr score of original feature
    d_state = mnrval(B{s},fet.data);
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


save(fullfile(Trial.spath,[Trial.filebase '-rec_lgr-randClassHeir-4.mat']), ...
              'iter_B',...
              'iter_bp_train',...
              'iter_bp_train_N',...
              'iter_bp_test',...
              'iter_bp_test_N',...
              'iter_new_smat')
ds = load(fullfile(Trial.spath,[Trial.filebase '-rec_lgr-randClassHeir-4.mat']));

ds

fileList = dir(Trial.spath);
files = {};
for file = fileList',
    if ~isempty(regexp(file.name,'rec_lgr')),
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


% hierarchy lies in the order of states
tnew_smat = new_smat;
for i=fliplr(1:size(new_smat,2)-1),
    tnew_smat(:,i+1) = new_smat(:,i+1).*double(~sum(new_smat(:,1:i),2));
end

LGRsmat = MTADxyz('data',tnew_smat,'sampleRate',fet.sampleRate);
smat = MTADxyz('data',~~stc2mat(Trial.stc,fet,states),'sampleRate',fet.sampleRate);


% cmat = zeros([numel(states),numel(states)]);
% for i = 1:numel(states),
%    for j = 1:numel(states),
%        cmat(i,j) = sum(double(smat(pper,i)==1&tnew_smat(pper,j)==1));
%    end
% end

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