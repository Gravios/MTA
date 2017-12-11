% MjgEd2016_supfig_FeatureCorregistrationTest
%
% Confirm estimation of feature offset is consistent across days given same rigidbody


trialList  ={'jg05-20120309.cof.all',...
             'jg05-20120310.cof.all',...
             'jg05-20120311.cof.all',...
             'jg05-20120312.cof.all',...
             'jg05-20120317.cof.all'
            };
nt = numel(trialList);

hfig = figure();
eds = linspace(-pi/2,pi/2,100);

for r = 1:numel(trialList),
    tempTrialList = trialList(1:nt~=r);    
    for t = 1:nt-1
        subplot2(nt,nt-1,r,t);  hold('on');
        Trial = MTATrial.validate(trialList{t});
        stc = Trial.load('stc','msnn_ppsvd_raux');
        pch = fet_HB_pitch(Trial);
        ind = stc{'x'};
        
        hax = bar(eds,histc(pch(ind,3),eds),'histc');
        hax.FaceAlpha = 0.5;
        hax.FaceColor = 'b';
        hax.EdgeColor = 'b';
        hax.EdgeAlpha = 0.5;
        
        map_to_reference_session(pch,Trial,trialList{r});
        
        hax = bar(eds,histc(pch(ind,3),eds),'histc');
        hax.FaceAlpha = 0.5;
        hax.FaceColor = 'r';
        hax.EdgeColor = 'r';
        hax.EdgeAlpha = 0.5;
        title({[Trial.filebase],['pitch ref: ' trialList{r}]});
    end
end

