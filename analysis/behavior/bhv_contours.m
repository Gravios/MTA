function hax = bhv_contours(varargin);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('sessionList',                'MjgER2016',                                      ...
                 'pitchReferenceTrial',        'Ed05-20140529.ont.all',                          ...
                 'states',                     {{'lloc+lpause&theta','hloc+hpause&theta',        ...
                                                 'rear&theta'}},                                 ...
                 'stateColors',                'wcr'                                             ...
);
[sessionList,pitchReferenceTrial,states,stateColors] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

errMsgs.badStateColors = 'MTA:analysis:behavior:bhv_contour:StateColorsMismatch';
errMsgs.badSessionList = 'MTA:analysis:behavior:bhv_contour:BadSessionList';


assert(numel(states)==numel(stateColors),errMsgs.badStateColors);


if ischar(sessionList),
    sessionList = get_session_list(sessionList);
elseif ~isstruct(sessionList),
    error(errMsgs.badSessionList);
end
numTrials = numel(sessionList);
numStates = numel(states);
pitchReferenceFeature = fet_HB_pitch(pitchReferenceTrial);
pitchReferenceTrial = repmat({pitchReferenceTrial},[1,numTrials]);
states = repmat({states},[1,numTrials]);



Trials = af(@(s)        MTATrial.validate(s),                 sessionList);
stc    = cf(@(t)        t.stc.copy(),                         Trials);
pch    = cf(@(t,rt,rf)  fet_HB_pitchB(t,[],[],[],rt,rf),  ...
            Trials,pitchReferenceTrial,repmat({pitchReferenceFeature},[1,numTrials]));

cf(@(p) set(p,'sampleRate',119.881035), pch );
cf(@(p) resample(p,10), pch );






for s = 1:numStates,
    hout{s} = cf(@(p,s,sts)  hist2(p(resample([s{sts}],p),:),linspace(-pi/2,pi/2,50),linspace(-pi/2,pi/2,50)), ...
                 pch,stc,repmat({states{s}},[1,numTrials]));
end

% DIAGNOSTIC 
figure,for t = 1:numel(hout{1}),clf();imagesc(hout{1}{t}');waitforbuttonpress();end

lpt = [stc{2}{'lloc+lpause&theta'}];
lpt.resample(pch{2});
figure,plot(pch{2}(lpt,1)),Lines(lpt.data(:),[],'k');

figure,
eds = linspace(-pi/2,pi/2,50);

edgs    = {eds};
edgs(2) = {eds};
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    


figure();
hax = gobjecs([1,3]);
for s = 1:numStates,
    hax(s) = subplot(1,numStates,s);
    sout = nansum(cat(3,hout{s}{:}),3);
    imagesc(eds,eds,sout');caxis([0,15000]);axis('xy');
    hold('on');
    o = conv2(sout,F,'same');
    contour(X,Y,o',[700,700],'linewidth',0.5,'Color',stateColors(s))
end
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]),  hax);



