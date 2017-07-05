


%Sessions = get_session_list('test_grp',...
%                '/storage/gravio/data/processed/xyz/',...
%                '/storage/gravio/data/processed/nlx/');

%% Model Info
train = true;
%states = {'walk','rear'};
states = {'walk','rear','sit','turn','shake','groom'};
fet = 'fet_lgr';
s.name = 'jg05-20120317';
s.trialName = 'all';
s.mazeName = 'cof';

model_names = {};
sWinner = 'r';
states =  states(~cellfun(@isempty,regexp(states,Trial.stc{sWinner}.label)));

for i = 1:numel(states),
%for s = Sessions
    Trial = MTATrial(s.name,s.trialName,s.mazeName);    
    model_names(end+1) = {[Trial.filebase,'-','pop_lgr-' states{i}]};%mfilename]};
    bhv_lgr(Trial,train,[states(i),Trial.stc{['a' sws]}.label],fet,model_names{end},false,true);
%end
end


ns = numel(Sessions);
pB = [];
for i = 1:numel(states),
%for s = 1:ns,
ds = load([model_names{i} '-fet_lgr-model.mat']);
pB = cat(3,pB,ds.B);
%end
end
figure,plot(d_state)
Lines(Trial.stc{'w'}(:),[],'g');
Lines(Trial.stc{'r'}(:),[],'r');

c = jet(ns);
figure,hold on,
for i = 1:ns,
    scatter(1:13,pB(:,2,i),30,c(i,:));
end

figure,imagesc(sq(pB(:,2,:)))