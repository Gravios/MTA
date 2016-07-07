function distributions_features(varargin)

% Default arguments for distributions_features(varargin)
defargs = struct(...
    'sList'          , 'hand_labeled'                                           ,...
    'featureSet'     , 'fet_mis'                                                ,...
    'featureOpts'    , {{'newSampleRate',12,'procOpts','SPLINE_SPINE_HEAD_EQD'}},...
    'featureName'    , ''                                                       ,...
    'state'          , 'gper'                                                   ,...
    'resample'       , true                                                     ,...
    'normalize'      , true                                                     ,...
    'mapToReference' , true                                                     ,...
    'RefTrial'       , 'jg05-20120317.cof.all'                                   ...
);


[sList,featureSet,featureOpts,featureName,state,...
 resample,normalize,mapToReference,RefTrial] = DefaultArgs(varargin,defargs,'--struct');

if isempty(featureName)
    featureName = featureSet;
end


sesList = SessionList(sList);

snames = cell([1,numel(sesList)]);

for s = 1:numel(sesList),
    snames{s} = sesList(s).sessionName; 
end

%Reference Trial Stuff

if mapToReference||normalize,
    RefTrial = MTATrial.validate(RefTrial);
end

if normalize,
    RefTrial = MTATrial.validate(RefTrial);
    RefState = RefTrial.stc{'a'};
    rfet = feval(featureSet,RefTrial,featureOpts{:});
    [rfet,Rmean,Rstd] = unity(rfet,[],[],[],[]);
end
    

cfet = {};
oStc = {};
for s = 1:numel(sesList),
    Trial = MTATrial.validate(sesList(s));
    Trial.load('stc',sesList(s).stcMode);
    if s == 1,
        [cfet{s},fett,fetd] = feval(featureSet,Trial,featureOpts{:});
    else
        [cfet{s}] = feval(featureSet,Trial,featureOpts{:});
    end

    if mapToReference&&~strcmp(Trial.filebase,RefTrial.filebase), 
        cfet{s}.map_to_reference_session(Trial,RefTrial); 
    end
    if normalize,
        cfet{s} = unity(cfet{s},[],Rmean,Rstd);           
    end

    Stc = Trial.load('stc',sesList(s).stcMode);

    if resample
        [Stc,~,cfet{s}] = resample_whole_state_bootstrap_noisy_trim(Stc,cfet{s},{state});
    end

    cfet{s}.data = cfet{s}(Stc{state},:);
end


if isempty(fett)||isempty(fetd)
    fett = repmat( {''}, [1,cfet.size(2)] );
    fetd = repmat( {''}, [1,cfet.size(2)] );
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
        hs.FaceAlpha = .3;
        cind = cind+1;
    end
    legend(snames{:});
    pause(.3)
    reportfig(fullfile(getenv('PROJECT'),'figures'), hfig, [featureName '-' state], 'features', false,sList, ...
              ['feature: ',num2str(f),' - ' fetd{f}],[],false,'png');
    close(hfig)           
end


