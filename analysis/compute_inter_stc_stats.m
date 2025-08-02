function labelingStats = compute_inter_stc_stats(sessionList,stcHL,stcCL,states,sampleRate);
%function labelingStats = compute_inter_stc_stats(sessionList,stc,states,sampleRate);
%
%

% if iscell(sessionList) && isa(sessionList{1},'MTAStateCollection'), blah blah,end


if isa(sessionList,'MTASession'),
% ENCAPSULATE Trial in cell
    Trials = {sessionList};
elseif ischar(sessionList),
% LOAD Trials    
    Trials = af(@(t)  MTATrial.validate(t), get_session_list(sessionList));
end

assert( iscell(Trials) & ~isempty(Trials) & isa(Trials{1},'MTASession'),...
        ['MTA:analysis:',mfilename,':TrialFormatExpectsCellArray'])


numStates = numel(states);
numTrials = numel(Trials);

% LOAD xyz as a reference for synchronization
xyz    = cf(@(Trial)  Trial.load('xyz')       , Trials);


if isempty(stcHL),
    stcHL = cf(@(Trial) Trial.load('stc'),Trials);
elseif isa(stcHL,'MTAStateCollection'),
    stcHL = {stcHL};
end
% SELECT only states to be compared
cf(@(stc,states) set(stc,'states',stc(states{:})),...
   stcHL,repmat({states},[1,numTrials]));
% CONVERT state collection into state matrix for inter labeler comparison
shl = cf(@(s,x,sts) MTADxyz('data',double(0<stc2mat(s,x,sts)),...
                            'sampleRate',x.sampleRate),...
         stcHL, xyz, repmat({states},[1,numTrials]));




if ischar(stcCL),
    stcCL = cf(@(t,s)  t.load('stc',s),  Trials,repmat({stcCL},[1,numTrials]));
elseif isa(stcCL,'MTAStateCollection')
    stcCL = {stcCL};
end
% SELECT only states to be compared
cf(@(stc,states) set(stc,'states',stc(states{:})),...
   stcCL,repmat({states},[1,numTrials]));
% CONVERT state collection into state matrix for inter labeler comparison
ysm = cf(@(s,x) MTADxyz('data',double(0<stc2mat(s,x)),'sampleRate',x.sampleRate),...
         stcCL,xyz);
for s = 1:numTrials,     ysm{s}.data(~nniz(xyz{s}),:) = 0;    end

% INITIALIZE labeling timeperiods
labelingEpochs = cf(@(t,x)  resample(t.stc{'a'}.cast('TimeSeries'),x), Trials,xyz);

% INITIALIZE inter labeler stats
labelingStats = repmat(struct('confusionMatrix',zeros([numStates,numStates]),...
                               'precision',     zeros([1,numStates]),...
                               'sensitivity',   zeros([1,numStates]),...
                               'accuracy',      0 ),[1,numTrials]);

% COMPUTE inter labeler stats if valid comparative Stc was loaded above
ind = cf(@(s,y,l) any(s.data,2)&any(y.data,2)&l.data==1, shl,ysm,labelingEpochs);
tcm = cf(@(s,y,i) confmat(s(i,:),y(i,:)),      shl,ysm,ind); % #DEP: netlab
af(@(l,t,s) setfield(l,'confusionMatrix',round(t{1}./s{1},2)),...
   labelingStats,tcm,repmat({sampleRate},[1,numTrials]));
for s = 1:numTrials,
    labelingStats(s).confusionMatrix = round(tcm{s}./sampleRate,2);
    labelingStats(s).precision = round(diag(tcm{s})'./sum(tcm{s},2)',4).*100;
    labelingStats(s).sensitivity = round(diag(tcm{s})'./sum(tcm{s}),4).*100;
    labelingStats(s).accuracy = sum(diag(tcm{s}))/sum(tcm{s}(:));
end
