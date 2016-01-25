function distributions_features(varargin)
[sList,featureSet,state,normalize,mapToReference,RefTrial] = DefaultArgs(varargin,{'hand_labeled_jg','fet_tsne_rev10','gper',true,false,{'jg05-20120317','all','cof'}});


sesList = SessionList(sList);

snames = cell([1,numel(sesList)]);

for s = 1:numel(sesList),
    snames{s} = sesList(s).sessionName; 
end



NEW_SAMPLE_RATE = 10;
RP_PATH = '/storage/gravio/figures';

%Reference Trial Stuff
if normalize,
    RefTrial = MTATrial(RefTrial{:});
    RefState = RefTrial.stc{'a'};
    rfet = feval(featureSet,RefTrial,NEW_SAMPLE_RATE);
    [rfet,Rmean,Rstd] = unity(rfet,[],[],[],[]);
end
    

cfet = {};
Stc = {};
for s = 1:numel(sesList),
    Trial = MTATrial(sesList(s).sessionName,...
                     sesList(s).trialName,...
                     sesList(s).mazeName);
    Trial.load('stc',sesList(s).stcMode);
    if s ==1,
        [cfet{s},fett,fetd] = feval(featureSet,Trial,NEW_SAMPLE_RATE);
    else
        [cfet{s}] = feval(featureSet,Trial,NEW_SAMPLE_RATE);
    end

    if mapToReference, cfet{s}.map_to_reference_session(Trial,RefTrial); end
    if normalize,      cfet{s} = unity(cfet{s},[],Rmean,Rstd);           end

    Stc = Trial.load('stc',sesList(s).stcMode);
    cfet{s}.data = cfet{s}(Stc{state},:);
end

c = jet(numel(cfet));
cind= 1;
for f = 1:cfet{1}.size(2),
    cind = 1;
    hfig = figure;
    hold('on');
    eds = linspace(-4,4,100);
    for s = 1:numel(cfet),
        hs = bar(eds,histc(cfet{s}(:,f),eds),'histc');
        hs.FaceColor = c(cind,:);    
        hs.FaceAlpha = .5;
        cind = cind+1;
    end
    legend(snames{:});
    pause(.3)
    reportfig(RP_PATH, hfig, featureSet, 'features', false,sList, ...
              ['feature: ',num2str(f)],[],false,'png');
    close(hfig)           
end


