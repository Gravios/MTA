


%Sessions = get_session_list('test_grp',...
%                '/storage/gravio/data/processed/xyz/',...
%                '/storage/gravio/data/processed/nlx/');

%% Model Info
train = true;
states = {'walk','rear','sit','turn','shake','groom'};
init_ns = numel(states);

Trial = MTATrial('jg05-20120317','all','cof');    
fet = 'fet_tsne';
model_names = {};
State_cat_order = {};%repmat({''},[init_ns,1]);
lgrm = {};
pB = repmat({[]},[init_ns,1]);


for s = 1:init_ns,
    disp(['inter: ' num2str(s) ', finding best state'])
    % 
    temp_states =  states(cellfun(@isempty,regexp(states,['(',strjoin(State_cat_order,'|'),')'])));
    
    if ~isempty(State_cat_order),
        sws = ['-' strjoin(State_cat_order,'-')];
    else
        sws = '';
    end
    
    for i = 1:numel(temp_states),
        model_names(s,i) = {[Trial.filebase,'-','pop_lgr-' temp_states{i} sws]};%mfilename]};
        bhv_lgr(Trial,...
                train,...
                [temp_states(i),Trial.stc{['a-' temp_states{i} sws]}.label],...
                fet,model_names{s,i},false,false);
    end

    
    for i = 1:numel(temp_states),
        lgrm{s,i} = load([model_names{s,i} '-fet_lgr-model.mat']);
        pB{s} = cat(2,pB{s},lgrm{s,i}.B);
    end
    
    [~,BestStateInd] = max(max(abs(pB{s}(2:end,:))));
    State_cat_order{s} = temp_states{BestStateInd};
    
end

figure,
imagesc(abs(pB{1}'));
colormap jet;




%% Periods known as Mobile
Trial = MTATrial('jg05-20120317');
Trial.load('stc','hand_labeled_rev2');
mper = Trial.stc{'w+n'};
mper.fillgaps;
mper.filename = 'jg05-20120317.cof.all.sst.mobile.b.mat';
mper.key = 'b';
mper.label = 'mobile';
mper.path = Trial.spath;
mper.save(1);

%lr fet god noooooooooooo
%% Test the by partition training
%QuickTrialSetup('jg05-20120317','lgr_training_set','cof',[20,0],[],1:3);
%QuickTrialSetup('jg05-20120317','lgr_testing_set','cof',[20,0],1:3);

Trial = MTATrial('jg05-20120317');
sampleRate = 30;
ofet = fet_lgr(Trial,sampleRate);
%ofet.resample(sampleRate);

%Trial = MTATrial('jg05-20120317','lgr_training_set');
%fet = fet_exp001(Trial);
fet = ofet.copy;
%fet = Trial.resync(fet);
%fet = fet_lgrN(Trial);




Trial.load('stc','hand_labeled_rev1');
aper = Trial.stc{'a'};
Trial.stc{'mobile'};
states = {'rear','mobile','groom','sit'};
B = {};
new_smat = [];
for s = states
    disp(['inter: ' s{1} ', finding best state'])

    % compute lgr coefficients
    [smat] = max(stc2mat(Trial.stc,fet,{s}),[],2);
    ind = resample(cast(aper.copy,'TimeSeries'),fet);
    ind = logical(ind.data);
    [B{end+1}] = mnrfit(fet(ind,:),smat(ind)+1,'model','nominal');

    % compute lgr score of original feature
    d_state = mnrval(B{end},fet.data);
    [~,dind] = max(abs(d_state),[],2);
    new_smat(:,end+1) = dind-1;
 
    % remove periods from aper which were scored as target behavior
    nper = ThreshCross(dind-1,.5,10);
    aper = aper-[Trial.stc{s{1}}];%nper
end


% hierarchy lies in the order of states
tnew_smat = new_smat;
for i=fliplr(1:size(new_smat,2)-1),
    tnew_smat(:,i+1) = new_smat(:,i+1).*double(~sum(new_smat(:,1:i),2));
end

%states = {'rear','mobile','groom','sit'};
smat = ~~stc2mat(Trial.stc,fet,states);

cmat = zeros([numel(states),numel(states)]);
for i = 1:numel(states),
    for j = 1:numel(states),
        cmat(i,j) = sum(double(smat(:,i)==1&tnew_smat(:,j)==1));
    end
end

cmat
bsxfun(@rdivide,diag(cmat),sum(cmat,1)')
bsxfun(@rdivide,diag(cmat),sum(cmat,2))

% 
% itr = 3;
% states_perm{itr} = states;
% B_perm{itr} = B;
% sacc_perm{itr} = sacc;
% cmat_perm{itr} = cmat;
% 
% 
% 
% Data.fillgaps(.2);


%T = MTATrial('jg05-20120310');
T = MTATrial('jg05-20120317','lgr_testing_set');
T.load('stc','hand_labeled_rev2');
T.stc{'turn+walk'};
tfet = ofet.copy;
tfet = T.resync(tfet);




test_smat = [];
for i = 1:numel(states),
    d_state = mnrval(B{i},tfet.data);
    [~,dind] = max(d_state,[],2);
    test_smat(:,end+1) = dind-1;
end
for i=fliplr(1:size(test_smat,2)-1),
    test_smat(:,i+1) = test_smat(:,i+1).*double(~sum(test_smat(:,1:i),2));
end

%figure,imagesc(test_smat')


tsmat = ~~stc2mat(T.stc,tfet,states);

tcmat = zeros([numel(states),numel(states)]);
for i = 1:numel(states),
    for j = 1:numel(states),
        tcmat(i,j) = sum(double(tsmat(:,i)==1&test_smat(:,j)==1));
    end
end

tsacc = bsxfun(@rdivide,diag(tcmat),sum(tcmat,2));


%% Test trained lgr model on a different session

T = MTATrial('jg05-20120310');
sfet = fet_lgr(T);


stest_smat = [];
for i = 1:numel(states),
    d_state = mnrval(B{i},sfet.data);
    [~,dind] = max(d_state,[],2);
    stest_smat(:,end+1) = dind-1;
end

for i=fliplr(1:size(stest_smat,2)-1),
    stest_smat(:,i+1) = stest_smat(:,i+1).*double(~sum(stest_smat(:,1:i),2));
end

states = {'rear','mobile','groom','sit'};
keys   = {   'r',     'b',    'm',  's'};
Stc = MTAStateCollection(T.spath,T.filebase,'lgr_rwms',T.stc.sync.copy,T.stc.origin);
Stc.states = {};
for i = 1:numel(states),
    per = ThreshCross(stest_smat(:,i),.5,5);
    per = bsxfun(@minus,per,[-1,0]);
    Stc.addState(Stc.path,T.filebase,per,sfet.sampleRate,Stc.sync.copy,sfet.origin,states{i},keys{i},'TimePeriods');
end
Stc.save(1)


stest_smat = tnew_smat;
states = {'rear','mobile','groom','sit'};
keys   = {   'r',     'b',    'm',  's'};
Stc = MTAStateCollection(Trial.spath,Trial.filebase,'lgr_rbms_rev3',Trial.stc.sync.copy,Trial.stc.origin);
Stc.states = {};
for i = 1:numel(states),
    per = ThreshCross(stest_smat(:,i),.5,5);
    per = bsxfun(@minus,per,[-1,0]);
    Stc.addState(Stc.path,Trial.filebase,per,fet.sampleRate,Stc.sync.copy,fet.origin,states{i},keys{i},'TimePeriods');
end
Stc.save(1)



%% Status - FAILED: LGR with mean values from bhv periods. 
% states = {'rear','turn+walk','groom','sit'};
% Bp = {};
% new_smat = [];
% for s = states
%     smv = [];
%     st = [Trial.stc{s{1},fet.sampleRate}];
%     for i = st.data',
%         smv(end+1,:) = mean(fet(i',:));
%     end
%     
%     ns = find(cellfun(@isempty,regexpi(states,s))==1);
%     
%     amv = [];
%     for i =1:numel(ns),
%         st = [Trial.stc{states{i},fet.sampleRate}];
%         for j = st.data'
%             amv(end+1,:) = mean(fet(j',:));
%         end
% 
%     end
% 
%     
p% % compute lgr coefficients
% smat = zeros([size(amv,1),1]);
% smat = cat(1,smat,ones([size(smv,1),1]));
% [Bp{end+1}] = mnrfit([amv;smv],smat+1,'model','nominal');
% 
% % compute lgr score of original feature
% d_state = mnrval(Bp{end},fet.data);
% [~,dind] = max(d_state,[],2);
% new_smat(:,end+1) = dind-1;
% 
% 
% % remove periods from aper which were scored as target behavior
% states(cellfun(@isempty,regexpi(states,s))==0)=[];
% %nper = ThreshCross(dind-1,.5,10);
% %aper = aper-nper;
% end
% 
% 
% tnew_smat = new_smat;
% for i=fliplr(1:size(new_smat,2)-1),
%     tnew_smat(:,i+1) = new_smat(:,i+1).*double(~sum(new_smat(:,1:i),2));
% end
% 
% smat = ~~stc2mat(Trial.stc,fet,states);
% 
% cmat = zeros([numel(states),numel(states)]);
% for i = 1:numel(states),
%     for j = 1:numel(states),
%         cmat(i,j) = sum(double(smat(:,i)==1&tnew_smat(:,j)==1));
%     end
% end
% 
% sacc = bsxfun(@rdivide,diag(cmat),sum(cmat,2));


%for i=1:4,Stc.states{i}.fillgaps;end

%states = {'rear','mobile','groom','sit'};
%[STC,d_state] = bhv_lgr(Trial,false,states,fet);

%LGRsmat = MTADxyz('data',~~stc2mat(STC,fet,states),'sampleRate',fet.sampleRate);
LGRsmat = MTADxyz('data',~~stc2mat(Stc,fet,states),'sampleRate',fet.sampleRate);
smat    = MTADxyz('data',~~stc2mat(Trial.stc,fet,states),'sampleRate',fet.sampleRate);
cmat = zeros([numel(states)+1,numel(states)+1]);
ind = Trial.stc{'a'};
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



[U,S,V] = svd(cov(fet(Trial.stc{'a-r-b'},:)));


figure,hist2(fet(Trial.stc{'m'},:)*V(:,1:2),100,100)

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2);
ang = create(MTADang,Trial,xyz);

figure,hold on
ind = Trial.stc{'a-r-b-m'};
edgs = linspace(0,pi,100);
bar(edgs,histc(abs(circ_dist(ang(ind,1,2,1),...
                                          ang(ind,4,7,1))),edgs),'histc');
ind = Trial.stc{'m'};
h = bar(edgs,histc(abs(circ_dist(ang(ind,1,2,1), ...
                                              ang(ind,4,7,1))),edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

figure
%ind = Trial.stc{'a-r-b'};
ind = Trial.stc{'a-r-b-m'};
ind = Trial.stc{'m'};
edgs ={};
edgs{1} = linspace(0,pi,100);
edgs{2} = linspace(.4,pi/2,100);
hist2([abs(circ_dist(ang(ind,1,2,1),ang(ind,4,7,1))),...
       ang(ind,1,2,2)],edgs{1},edgs{2})
caxis([0,100])


figure,hold on
ind = Trial.stc{'a-r-b-m'};
bar(linspace(-pi*2,pi*2,100),...
histc(...
circ_dist(ang(ind,1,2,1),ang(ind,2,3,1))+...
circ_dist(ang(ind,2,3,1),ang(ind,3,4,1))+...
circ_dist(ang(ind,4,5,1),ang(ind,5,7,1)),...
linspace(-pi,pi*2,100)),'histc');
ind = Trial.stc{'m'};
h = bar(linspace(-pi*2,pi*2,100),...
histc(...
circ_dist(ang(ind,1,2,1),ang(ind,2,3,1))+...
circ_dist(ang(ind,2,3,1),ang(ind,3,4,1))+...
circ_dist(ang(ind,4,5,1),ang(ind,5,7,1)),...
linspace(-pi,pi*2,100)),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;


figure,hold on
ind = Trial.stc{'a-r-m'};
edgs = linspace(-6,1,100);
bar(edgs,...
histc(...
log10(abs(mean([circ_dist(ang(ind,1,2,1),ang(ind,2,3,1)),...
circ_dist(ang(ind,1,3,1),ang(ind,3,4,1)),...
circ_dist(ang(ind,1,5,1),ang(ind,5,7,1))],2))),...
edgs),'histc');
ind = Trial.stc{'m'};
h = bar(edgs,...
histc(...
log10(abs(mean([circ_dist(ang(ind,1,2,1),ang(ind,2,3,1)),...
circ_dist(ang(ind,1,3,1),ang(ind,3,4,1)),...
circ_dist(ang(ind,1,5,1),ang(ind,5,7,1))],2))),...
edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;



figure,hold on
ind = Trial.stc{'a-r-m'};
edgs = linspace(-2,1,100);
bar(edgs,...
histc(...
log10(abs(mean([circ_dist(ang(ind,1,2,2),ang(ind,2,3,2)),...
circ_dist(ang(ind,1,3,2),ang(ind,3,4,2)),...
circ_dist(ang(ind,1,5,2),ang(ind,5,7,2))],2))),...
edgs),'histc');
ind = Trial.stc{'m'};
h = bar(edgs,...
histc(...
log10(abs(mean([circ_dist(ang(ind,1,2,2),ang(ind,2,3,2)),...
circ_dist(ang(ind,1,3,2),ang(ind,3,4,2)),...
circ_dist(ang(ind,1,5,2),ang(ind,5,7,2))],2))),...
edgs),'histc');
h.FaceColor = 'r';
h.FaceAlpha = .5;

scur = log10(abs(mean([circ_dist(ang(:,1,2,1),ang(:,2,3,1)),...
                       circ_dist(ang(:,1,3,1),ang(:,3,4,1)),...
                       circ_dist(ang(:,1,5,1),ang(:,5,7,1))],2)));

scurp = log10(abs(mean([circ_dist(ang(:,1,2,2),ang(:,2,3,2)),...
                       circ_dist(ang(:,1,3,2),ang(:,3,4,2)),...
                       circ_dist(ang(:,1,5,2),ang(:,5,7,2))],2)));

scur = MTADxyz('data',scur,'sampleRate',ang.sampleRate);
scurp = MTADxyz('data',scurp,'sampleRate',ang.sampleRate);


ind = Trial.stc{'a-r-b-m'};

ind = Trial.stc{'m'};
edgs = {};
edgs{1} = linspace(-2,.5,100);
edgs{2} = linspace(-.8,.2,100);
figure
hist2([scur(ind),scurp(ind)],edgs{1},edgs{2}),
caxis([0,400])


ind = ':';%Trial.stc{'a'}
figure,plot(mean([circ_dist(ang(ind,1,2,1),ang(ind,2,3,1)),...
circ_dist(ang(ind,1,3,1),ang(ind,3,4,1)),...
circ_dist(ang(ind,1,5,1),ang(ind,5,7,1))],2))

Lines(Trial.stc{'m'}(:),[],'g');

%% Validation

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
    bp_train{s}.cast('TimeSeries');
    bp_train{s}.resample(fet);

    % compute lgr coefficients, 
    %atemper = atemper-bp_test{s};
      %bypass for parfor . . . the stupid @!#$
      bp_test{s}.resample(atemper.sampleRate);
      atemper.data = SubstractRanges(atemper.data,bp_test{s}.data);
      atemper.label = [atemper.label '-' bp_test{s}.label];
      atemper.key = '';

    
    ind = resample(cast(atemper.copy,'TimeSeries'),fet);
    ind = logical(ind.data);

    [B{s}] = mnrfit(fet(ind,:),bp_train{s}(ind)+1,'model','nominal');
    % compute lgr score of original feature
    d_state = mnrval(B{s},fet.data);
    [~,dind] = max(abs(d_state),[],2);
    new_smat(:,s) = dind-1;
    % remove periods from aper which were scored as target behavior
    nper = ThreshCross(dind-1,.5,10);
    %aper = aper-sper;
      %bypass for parfor . . . the stupid @!#$, every !@#$#@ minus sign    
      sper.resample(aper.sampleRate);
      aper.data = SubstractRanges(aper.data,sper.data);
      aper.label = [aper.label '-' sper.label];
      aper.key = '';

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

iters = 5;

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

nanmean(iter_ppv)
figure,hist(iter_ppv(:,1),edgs)
figure,hist(iter_ppv(:,2),edgs)
figure,hist(iter_ppv(:,3),edgs)
figure,hist(iter_ppv(:,4),


