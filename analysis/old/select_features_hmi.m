function [stateOrd,fetInds,miAstates] = select_features_hmi(Trial,stc,fet,varargin)
% function [stateOrd,fetInds] = select_features_hmi(Trial,stc,fet)
% [states,display] = DefaultArgs(varargin,{{'rear','walk','turn','pause','groom','sit'},false});
% select best features for binary classifier of state vs all/state
% perform heirarchical analysis of mutual information gain at each
% level for each state vs all/state.
%
% Use d' to find best features to discriminate between the final 
% pair of states.

[states,display] = DefaultArgs(varargin,{{'rear','walk','turn','pause','groom','sit'},false});

fetInds = {};
stateOrd = {};
miAstates = {};
gStates = states;
fetRanges = {};

while numel(gStates) > 2

    [mixy,fetRanges] = calculate_MI_states_vs_features(stc,fet,gStates,fetRanges);

% ACCUMULATE the positive changes in mutual information between all states and minus one state.
    dms = [];
    for s = 1:numel(gStates)
        dm = (mixy(1,:)-mixy(s+1,:))';
        scatter(s,sum(dm(dm>0)),20);
        dms(end+1) = sum(dm(dm>0));
    end

% FIND the state which with the greatest positive chang
    [~,sind] = max(dms);
    dm = (mixy(1,:)-mixy(sind+1,:))';
% FIND the feature indices where the change in mutual information is greater than 0.2bits;
    fetInds{end+1} =  find(dm>0.20);
    stateOrd{end+1} = gStates{sind};

    if display,
        hfig = figure(283823899);
        eds = linspace(-.25,.75,1000);
        for s = 1:numel(gStates),
            subplot(2,ceil(numel(gStates)/2),s);
            plot(mixy([1,s+1],:)');
            bar(eds,histc(mixy(1,:)-mixy(s+1,:),eds),'histc');
            title({'Distribution of states-feature MI difference between',...
                  ['all states VS ' gStates{s} ' held out']});
            ylabel('Count');
            xlabel('Difference of Mutual Information (bits)');
        end

        reportfig(fullfile(Trial.path.project,'figures'),... 
                  hfig,...          Figure Handle 
                  strjoin({Trial.filebase,fet.label,stc.mode,gStates{:}},'-'),... FileName
                  mfilename,...     Subdirectory 
                  false,...         Preview
                  strjoin({Trial.filebase,fet.label,stc.mode,gStates{:}},'-'),... Tag
                  strjoin({Trial.filebase,fet.label,stc.mode,gStates{:}},'-'),... Comment
                  200,...           Resolution
                  true,...          SaveFig
                  'png',...         Format
                  16,10);%          width & height (cm)
    end
    miAstates{end+1} = mixy;
    gStates(~cellfun(@isempty,regexp(gStates,['^',stateOrd{end},'$'])))=[];

end

% Use d-prime to find best features to separate final two features
dp = (nanmean(fet(stc{gStates{1}},:))-nanmean(fet(stc{gStates{2}},:)))./...
     (0.5*(nanvar(fet(stc{gStates{1}},:))+nanvar(fet(stc{gStates{2}},:))));

fetInds{end+1} = find(abs(dp)>2)';
stateOrd = cat(2,stateOrd,gStates); % Need to specify order at some point


