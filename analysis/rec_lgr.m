


%Sessions = SessionList('test_grp',...
%                '/storage/gravio/data/processed/xyz/',...
%                '/storage/gravio/data/processed/nlx/');

%% Model Info
train = true;
states = {'walk','rear','sit','turn','shake','groom'};
init_ns = numel(states);

Trial = MTATrial('jg05-20120317','all','cof');    
fet = fet_lgrN(Trial);
fet.resample(30);
model_names = {};
State_cat_order = {};%repmat({''},[init_ns,1]);
lgrm = {};
pB = repmat({[]},[init_ns,1]);

 
for s = 1:init_ns,
    disp(['inter: ' num2str(s) ', finding best state'])
    temp_states =  states(cellfun(@isempty,regexp(states,['(',strjoin(State_cat_order,'|'),')'])));
    
    if ~isempty(State_cat_order),
        sws = ['-' strjoin(State_cat_order,'-')];
    else
        sws = '';
    end
    
    for i = 1:numel(temp_states),
        model_names(s,i) = {[Trial.filebase,'-','pop_lgr-' temp_states{i} sws]};%mfilename]};
        bhv_lgr(Trial,train,[temp_states(i),Trial.stc{['a-' temp_states{i} sws]}.label],...
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






%lr fet god noooooooooooo
%% Test the by partition training
%QuickTrialSetup('jg05-20120317','lgr_training_set','cof',[20,0],[],1:3);
%QuickTrialSetup('jg05-20120317','lgr_testing_set','cof',[20,0],1:3);

Trial = MTATrial('jg05-20120317');
sampleRate = 30;
ofet = fet_lgr(Trial);
ofet.resample(sampleRate);

%Trial = MTATrial('jg05-20120317','lgr_training_set');
%fet = fet_exp001(Trial);
fet = ofet.copy;
fet = Trial.resync(fet);
%fet = fet_lgrN(Trial);


Trial.load('stc','hand_labeled_rev2');
aper = Trial.stc{'a'};
% Initialize state since stc2mat has a bug
Trial.stc{'turn+walk'};


states = {'rear','turn+walk','groom','sit'};
B = {};
new_smat = [];
mask = false([fet.size(1),1]);
mask(randi(
for s = states
        disp(['inter: ' s{1} ', finding best state'])

% compute lgr coefficients
[smat] = max(stc2mat(Trial.stc,fet,{s}),[],2);
ind = resample(cast(aper.copy,'TimeSeries'),fet);
ind = logical(ind.data);
[B{end+1}] = mnrfit(fet(ind,:),smat(ind)+1,'model','nominal');

% compute lgr score of original feature
d_state = mnrval(B{end},fet.data);
[~,dind] = max(d_state,[],2);
new_smat(:,end+1) = dind-1;

% remove periods from aper which were scored as target behavior
nper = ThreshCross(dind-1,.5,10);
aper = aper-nper;
end
% 
% d_sit = mnrval(B{4},fet.data);
% B_sit = B{4}(2:end);
% f_sit = 1;%find(abs(B_sit)>1);
% 
% B_sel = {};
% [smat] = max(stc2mat(Trial.stc,fet,{'walk'}),[],2);
% smat = ~~smat;
% ind = resample(cast(aper.copy,'TimeSeries'),fet);
% ind.data(isnan(ind.data)) = 0;
% ind = logical(ind.data);
% [B_sel{end+1}] = mnrfit(fet(ind,f_sit),smat(ind)+1,'model','nominal');
% 
% d_sit_sel = mnrval(B_sel{1},fet.data(:,f_sit));
% 
% figure,hold on,plot(d_sit),plot(d_sit_sel)

% hierarchy lies in the order of states
tnew_smat = new_smat;
for i=fliplr(1:size(new_smat,2)-1),
    tnew_smat(:,i+1) = new_smat(:,i+1).*double(~sum(new_smat(:,1:i),2));
end

smat = ~~stc2mat(Trial.stc,fet,states);

cmat = zeros([numel(states),numel(states)]);
for i = 1:numel(states),
    for j = 1:numel(states),
        cmat(i,j) = sum(double(smat(:,i)==1&tnew_smat(:,j)==1));
    end
end

sacc = bsxfun(@rdivide,diag(cmat),sum(cmat,2));
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
% % compute lgr coefficients
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





