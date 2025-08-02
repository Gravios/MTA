function distributions_features(varargin)
% Default arguments for distributions_features(varargin)
% defargs = struct(...
%     'sList'          , 'hand_labeled'                                           ,...
%     'featureSet'     , 'fet_mis'                                                ,...
%     'featureOpts'    , {{'newSampleRate',12,'procOpts','SPLINE_SPINE_HEAD_EQD'}},...
%     'featureName'    , ''                                                       ,...
%     'state'          , 'gper'                                                   ,...
%     'resample'       , true                                                     ,...
%     'normalize'      , true                                                     ,...
%     'mapToReference' , true                                                     ,...
%     'RefTrial'       , 'jg05-20120317.cof.all'                                   ...
% );



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(                                                                                ...
    'sessionList',              'hand_labeled',                                                  ...
    'featureSet',                'fet_mis',                                                      ...
    'featureOpts',               {{'newSampleRate',12}},                                         ...
    'state',                     'pause',                                                         ...
    'randomizationMethod',       @resample_whole_state_bootstrap_noisy_trim,                     ...
    'normalize',                 true,                                                           ...
    'map2reference',             true,                                                           ...
    'normalizationSessionList',  'hand_labeled',                                                 ...    
    'referenceTrial',            'jg05-20120317.cof.all'                                         ...
);

[sessionList,featureSet,featureOpts,state, randomizationMethod,                                  ...
normalize,map2reference,normalizationSessionList,referenceTrial] =                               ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% LOAD Trials 
Trials = af(@(t)  MTATrial.validate(t),  get_session_list(sessionList));
snames = cf(@(t)  t.filebase,            Trials);
numTrials = numel(Trials);

% LOAD features 
[features,fett,fetd] = cf(@(Trial,fetSet,fetOpts) feval(fetSet,Trial,fetOpts{:}),...
                          Trials,...
                          repmat({featureSet },[1,numTrials]),...
                          repmat({featureOpts},[1,numTrials]));

% MAP features to reference session
if map2reference,
    xyzls  = cf(@(Trial)  Trial.load('xyz'),        Trials);
             cf(@(f,x)    x.resample(f),            features,xyzls);
    refTrial = MTATrial.validate(referenceTrial);
    cf(@(f,t,r) f.map_to_reference_session(t,r),...
       features, Trials, repmat({refTrial},[1,numTrials]));
    for s = 1:numTrials, features{s}.data(~nniz(xyzls{s}),:,:) = 0;end

    if normalize,
% CONDITIONAL normalization, use multiple sessions for normalization        
        [refMean,refStd] = load_normalization_parameters_unity(featureSet,...
                                                               refTrial.filebase,...
                                                               normalizationSessionList);
        cf(@(f,m,s) f.unity(@nan,m,s), ...
           features,...
           repmat({refMean},[1,numTrials]),...
           repmat({refStd}, [1,numTrials]));
    end

elseif normalize,
% CONDITIONAL normalization, use current sessions for normalization    
    zfrCat = cf(@(f) get(f,'data'),    features); 
    zfrCat = cat(1,zfrCat{:});
    zfrMean = nanmean(zfrCat(nniz(zfrCat),:,:));
    zfrStd  = nanstd( zfrCat(nniz(zfrCat),:,:));
    cf(@(w,m,s) set(w,'data',nunity(w,[],m,s)),...
       features,...
       repmat({zfrMean},[1,numTrials]),...
       repmat({zfrStd}, [1,numTrials]));
    clear('zfrCat','zfrMean','zfrStd');
end


% LOAD state collection
Stc = cf(@(t) t.stc.copy(), Trials);

if ~isempty(randomizationMethod),
% RESAMPLE randomly the features from the states
    [Stc,~,features] = cf(@(s,f,sts) randomizationMethod(s,f,sts), ...
                          Stc,features,repmat({{'walk','turn','pause'}},[1,numTrials]));
end

cf(@(c,s,sts)  set(c,'data',c(s{sts},:)),  features, Stc, repmat({state},[1,numTrials]));


% SET feature labels and descriptions
if isempty(fett)||isempty(fetd)
    fett = repmat( {''}, [1,size(features{1},2)] );
    fetd = repmat( {''}, [1,size(features{1},2)] );
else
    fett = fett{1};
    fetd = fetd{1};
end


% PLOT feature distributions
c = jet(numTrials);
for f = 1:size(features{1},2),
    hfig = figure();
    hold('on');
    eds = linspace(-pi,pi,100);
    for s = 1:numTrials,
        subplot(numTrials,1,s);
        bar(eds,histc(features{s}(:,f),eds),'histc');
        xlim([-pi,pi]);
    end
    pause(.1)
    reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path
              hfig,                                   ... Figure Handel
              [features{1}.label,'-',                 ... Figure Name
               'N',num2str(normalize),                ... ...
               'M',num2str(map2reference),'-',state], ... ...
              'features',                             ... Figure Directory
              false,                                  ... Preview
              sessionList,                            ... Tag
              ['feature: ',num2str(f),' - ' fetd{f}], ... Comment
              [],                                     ... Resolution
              false,                                  ... Save Figure
              'png',                                  ... File Format
              8,                                      ... Width
              4,                                      ... Heigth
              f                                       ... Index
    );
end

% END MAIN -----------------------------------------------------------------------------------------