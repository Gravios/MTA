function [ixy,metadata] = compute_cross_correlation(varargin)

% TODO : featureSet Should accept an MTADfet object
% TODO : circle shifting can be replaced with subtraction based
%        index shifting, pulling segments at state index points

% SET the trial set examined in the anaylsis
% SET the feature set examined in the anaylsis
% SET subset of features tobe included in the analysis
% SET edges to domain of the log10 marker speed
% SET the behavioral state subset
% SET shift values from -2 seconds to 2 seconds

global MTA_PROJECT_PATH

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',              'hand_labeled',                                     ...
                 'featureSet',               'fet_bref',                                         ...
                 'featureSubset',            [16:2:24,17:2:25,26:30],                            ...
                 'state',                    'walk',                                             ...
                 'edges',                    linspace(-1,1,64),                                  ...
                 'sbound',                   -240:240,                                           ...
                 'sampleRate',               119.991035,                                         ...
                 'overwrite',                false,                                              ...
                 'tag',                      ''                                                  ...
);
[sessionList,featureSet,featureSubset,state,edges,sbound,sampleRate,                             ...
 overwrite,tag] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------

% LOAD Trials.
Trials = af(@(Trial) MTATrial.validate(Trial), get_session_list(sessionList));
% LOAD state collections.
stc = cf(@(t) t.load('stc'), Trials);

metadata = struct('sessionList',            sessionList,...
                  'featureSet',             featureSet,...
                  'featureSubset',          featureSubset,...
                  'stcMode',                {cf(@(s) s.mode, stc)},...
                  'state',                  {state},...
                  'edges',                  edges,...
                  'sbound',                 sbound,...
                  'sampleRate',             sampleRate);

% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash(metadata);
end
%---------------------------------------------------------------------------------------------------

% SET file path where analysis results are tobe saved
FileName = fullfile(MTA_PROJECT_PATH,'analysis',['xcorr-',tag,'.mat']);

if ~exist(FileName,'file') || overwrite,
% CONVERT speed to log10 scale
    features = cf(@(t) fet_bref(t), Trials);


% MAP to a reference session
    cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'), ...
       features, Trials);

% LOAD multi-session normalization parameters
    [refMean,refStd] = load_normalization_parameters_unity('fet_bref',...
                                                      'jg05-20120317.cof.all',...
                                                      'hand_labeled');
% NORMALIZE feature set
    cf(@(f,m,s) f.unity(@nan,m,s), ...
       features,...
       repmat({refMean},[1,numel(Trials)]),...
       repmat({refStd}, [1,numel(Trials)]));
    
% FILTER Feature set 
    cf(@(f)  f.filter('ButFilter',3,10,'low'),  features);

% LOAD multi-session feature domain estimates 
% $$$     psa = load_normalization_parameters_mapminmax(features{1}.label,...
% $$$                                                   'jg05-20120317.cof.all',... 
% $$$                                                   'hand_labeled');
% $$$     cf(@(f,p)    set(f,'data',mapminmax('apply',f.data',p)'), features,repmat({psa},[1,numel(Trials)]));

% REMOVE non-dynamic features. NOTE: lazy to remove after above transformations
    cf(@(f,sub)  set(f,'data',f(:,sub)),   features,repmat({featureSubset},[1,numel(Trials)]));

% UNENCAPSULATE data from MTA object
    v = cf(@(f)  f.data,  features);
    v = cat(1,v{:}); 

% INITIALIZE mutual information matrix
    ixy = zeros([numel(sbound),size(v,2),size(v,2)]);

    if ischar(state),
% SELECT periods of behavior and cast to TimeSeries
        ind = cf(@(s,sts)  cast([s{sts}],'TimeSeries'),  stc,repmat({state},[1,numel(Trials)]));
              cf(@(i,f)    resample(i,f)              ,  ind,features);
        ind = cf(@(i)      i.data                     ,  ind);
    
    elseif iscell(state),
% SELECT state spcific events
% SELECT time points of a specific transition type and cast to TimeSeries
        transIndex = cf(@(s,t,sts,f) get_state_transitions(s,t,sts,[],f), ...
                        stc,Trials,repmat({state},[1,numel(Trials)]),...
                        features);
        transIndex = cf(@(i)       round(mean(i,2)) ,  transIndex);
        %transIndex = cf(@(i) subsasgn(i,substruct('()',{find(diff([0;i]) < 120),{[]}})), transIndex);
        %for s = 1:numel(Trials), transIndex{s}(diff([0;i]) < 120) = []; end
        ind        = cf(@(f)       zeros([size(f,1),1]),  features);
        ind        = cf(@(i,trans) subsasgn(i,substruct('()',{trans,1}),1), ...
                        ind,transIndex);
% $$$         ind        = cf(@(i)       conv(i,ones([60,1]),'same')>=1,ind);
          
    elseif isa(state,'MTADepoch'),
        
    end

% SELECT random time points from a state subset
    
    ind = find(logical(cat(1,ind{:})));
    nind = sum(ind==1);
    sbound = [0,sbound];

    
    for m = 1:size(v,2)
        for o = 1:size(v,2)
            sv = v;
            for s = 1:numel(sbound)-1,
                sv = circshift(sv,-diff(sbound([s,s+1])));
                segs = prod(cat(3,GetSegs(v(:,m),ind,60,1),GetSegs(sv(:,o),ind,60,1)),3);
                segs(:,~nniz(segs')) = [];
                ixy(s,m,o) = sum(median(segs,2));
            end
        end
    end

    save(FileName,'ixy','metadata');
    return
end

load(FileName);
