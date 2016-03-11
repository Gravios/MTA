function [stateOrd,fetInds] = select_features_hmi(Trial,stc,fet)
[states,display] = DefaultArgs(varargin,{{'rear','walk','turn','pause','groom','sit'},false});

fetInds = {};
stateOrd = {};

while numel(gStates) > 2

    mixy = calculate_MI_states_vs_features(stc,fet,states);

    dms = [];
    for s = 1:numel(gStates)
        dm = (mixy(1,:)-mixy(s+1,:))';
        scatter(s,sum(dm(dm>0)),20);
        dms(end+1) = sum(dm(dm>0));
    end

    [~,sind] = max(dms);
    dm = (mixy(1,:)-mixy(sind+1,:))';
    fetInds{end+1} =  find(dm>0.20);
    staetOrd{end+1} = gStates{sind};

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

        reportfig(fullfile(Trial.path.data,'figures'),... 
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

    gStates(sind) = [];
end

