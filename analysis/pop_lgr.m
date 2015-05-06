


Sessions = SessionList('test_grp',...
                '/storage/gravio/data/processed/xyz/',...
                '/storage/gravio/data/processed/nlx/');

%% Model Info
train = true
states = {'walk','rear'};
fet = 'fet_lgr';

model_names = {};
for s = Sessions
    Trial = MTATrial(s.name,s.trialName,s.mazeName);    
    model_names(end+1) = {[Trial.filebase,'-','pop_lgr']};%mfilename]};
    bhv_lgr(Trial,train,states,fet,model_names{end},false,true);
end


ds = load([model_names{3} '-fet_lgr-model.mat'])
ds.B

figure,plot(d_state)
Lines(Trial.stc{'w'}(:),[],'g');
Lines(Trial.stc{'r'}(:),[],'r');