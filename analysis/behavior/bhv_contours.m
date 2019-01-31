function [C,H] = bhv_contours(varargin);
%function [C,hax] = bhv_contours(varargin);
% create contour map of feature jpdf
%
%
% only works for fet_HB_pitchB for the moment

global MTA_PROJECT_PATH

diagnostic = false;

% varargin = {};

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionListName',            'MjgER2016',                                      ...
                 'featureSet',                 'fet_HB_pitchB',                                  ...
                 'featureInd',                 [1,2],                                            ...
                 'featureBin',                 {{linspace(-2,2,50),linspace(-2,2,50)}},          ...
                 'referenceTrial',             'Ed05-20140529.ont.all',                          ...
                 'states',                     {{'lloc+lpause&theta','hloc+hpause&theta',        ...
                                                 'rear&theta'}},                                 ...
                 'stateColors',                'wcr',                                            ...
                 'tag',                        '',                                               ...
                 'overwrite',                  false                                             ...
);
[sessionListName,featureSet,featureInd,featureBin,referenceTrial,states,stateColors,tag,overwrite] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

errMsgs.badStateColors = 'MTA:analysis:behavior:bhv_contour:StateColorsMismatch';
errMsgs.badSessionList = 'MTA:analysis:behavior:bhv_contour:BadSessionList';

assert(numel(states)==numel(stateColors),...
       errMsgs.badStateColors,...
       'Number of colors does not match number of states');

try,
    sessionList = get_session_list(sessionListName);
catch err, disp(err);
    error(errMsgs.badSessionList,'The session list name provided was not found');
end



% TAG creation -------------------------------------------------------------------------------------
if isempty(tag),
    tag = DataHash({sessionList,...
                    featureSet,...
                    featureInd,...
                    featureBin,...                    
                    referenceTrial,...
                    states});
end
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

fileName = fullfile(MTA_PROJECT_PATH,'analysis',['bhv_contours','-',tag,'.mat']);
if ~exist(fileName,'file') || overwrite,
    
    numTrials = numel(sessionList);
    numStates = numel(states);

% LOAD Trials
% LOAD Behavior State Collections
% LOAD Feature Sets
    Trials = af(@(s)  MTATrial.validate(s),  sessionList);
    Stc    = cf(@(T)  T.stc.copy(),          Trials);

    
% REMOVE third session: error in stc--------
    if strcmp (sessionListName, 'MjgER2016')
        Trials(3:4) = [];
        Stc(3:4) = [];
    end
% REMOVE -----------------------------------


% LOAD features
    if ischar(featureSet) || ~all(isa(featureSet,'MTADfet')),
        Fet    = cf(@(T,f,r)  feval(f,T,[],[],[],r),  ...
                    Trials,...
                    repmat({featureSet},size(Trials)),...
                    repmat({referenceTrial},size(Trials)));
    else
        Fet = featureSet;
    end
    
    cf(@(p) set(p,'sampleRate',119.881035), Fet);

    
% ACCUMULATE 2d histograms
% NOTE may break due to the nonzero part of the anon func, no comment -_-
    for s = 1:numStates,
        hout{s} = cf(@(fetSet,stc,sts,fetInds,bins)                                        ...
                     hist2(reshape(nonzeros(fetSet(stc{[sts,'&gper']},fetInds)),[],2),      ...
                           bins{:}),                                                       ...
                     Fet,                                                                  ...
                     Stc,                                                                  ...
                     repmat({states{s}}, size(Trials)),                                    ...
                     repmat({featureInd},size(Trials)),                                    ...
                     repmat({featureBin},size(Trials)));
    end


% DIAGNOSTIC 
    if diagnostic,
        figure();  
        for t = 1:numel(hout{2}),  
            clf();  
            for s = 1:numStates,
                subplot(1,numStates,s);  imagesc(hout{s}{t}');  title(states{s})
                axis('xy');
            end            
            title(num2str(t))
            waitforbuttonpress();  
        end
        %for s = 1:3,hout{s}(3) = [];end
    end
    

% SET edges and smoothing parameters
    figure,
    edx = featureBin{1};
    edy = featureBin{2};
    edgs    = {edx};
    edgs(2) = {edy};
    [edgs{:}] = get_histBinCenters(edgs);
    [X,Y] = meshgrid(edgs{:});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];

% PLOT contours for each state
    hax = gobjects([1,3]);
    C = cell([1,3]);
    H = cell([1,3]);
    for s = 1:numStates,
        hax(s) = subplot(1,numStates,s);
        sout = nansum(cat(3,hout{s}{:}),3);
        imagesc(edx,edy,sout');caxis([0,2000]);axis('xy');
        hold('on');
        o = conv2(sout,F,'same');
        o = o/sum(o(:));
        oPrct = 0;
        oThresh = 0.1;
        wCounter = 0;
        while oPrct<0.8,
            oPrct = sum(o(o(:)>oThresh));
            oThresh = oThresh/1.1;
            if wCounter > 1200, break, end
        end
        [C{s},H{s}] = contour(X,Y,o',[oThresh,oThresh],'linewidth',0.5,'Color',stateColors(s));
    end

% SAVE contour objects
    save(fileName,'H','C');
else,
% LOAD contour objects
    load(fileName);
end

% END MAIN -----------------------------------------------------------------------------------------