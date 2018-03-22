function Stc = label_bhv_reduced_aux(Stc,Trial,varargin)
% function Trial = label_bhv_reduced_aux(Trial,varargin)
%
% Subclassification of locomotive and paused states into high and low sub-behaviors
%
% Stc:   (MTAStateCollection)   Behavior state collection 
%        (string)               contains the mode label for loading a collection
%
% Trial: (MTATrial)             Trial of Stc for labeling
%        (string)               validation reference in form (Session.exp.trial)
%                               (eg 'jg05-20120317.cof.all')
%
% varargin:
%    pitchTreshold             (Numeric)                    -0.65
%    referenceTrial            (String)                     'Ed05-20140529.ont.all'
%    overwrite                 (Logical)                    false 
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('pitchThreshold',         -0.65,                                                ...
                 'referenceTrial',         'Ed05-20140529.ont.all',                              ...
                 'overwrite',              false                                                 ...
);
[pitchThreshold,referenceTrial,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

if isempty(Stc), 
    Stc = Trial.stc.copy;
elseif ischar(Stc),
    Stc = Trial.load('stc',Stc);
end


% LOAD pitches 
% MAP pitches to reference trial
pch = fet_HB_pitch(Trial);
map_to_reference_session(pch,Trial,referenceTrial);    
filter(pch,'ButFilter',3,1.5,'low');


% SPLIT locomotion (loc.x) state into high and low
sempty = isempty(Stc{'h'})||isempty(Stc{'l'});
if sempty||overwrite,
% REMOVE existing states with same keys
    if ~sempty, 
        Stc.states(Stc.gsi('h')) = [];
        Stc.states(Stc.gsi('l')) = [];
    end
    


    if isempty(pitchThreshold)
% RETRIVE threshold from on pitch-rhm analysis
%         (MTA/analysis/rhythmic_head_motion/compute_session_rhm_distribution.m)
        fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
        if ~exist(fpath,'file')||overwrite,  
            Trial.stc = Stc.copy;
            compute_session_rhm_distribution(Trial,[],'loc','NN0317R'); 
        end
        afig = hgload(fpath);
        rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','compute_session_rhm_distribution-hangle'));
        rhm_distrb = rhm_distrb(1);
        figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
        mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
        mrhmp(mrhmp==0) = nan;
        rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
            +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
        pitchThreshold = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
        delete(afig);
    end    
    
% INITIALIZE new states' objects with copy
    locPeriods = Stc{'x'};
    wind = locPeriods.copy();
    lang = locPeriods.copy();
    hang = locPeriods.copy();    

% CREATE low loc state
    lang.data = ThreshCross(pch(:,3)<pitchThreshold,.5,1);
    lang.label = 'lloc';
    lang.key = 'l';
    lang.sampleRate = pch.sampleRate;
    lang.filename = [Trial.filebase '.sst.lloc.l.mat'];
    lang.resample(locPeriods.sampleRate);
    Stc.states{end+1} = lang&wind.data;


% CREATE high loc state
    hang.data = ThreshCross(pch(:,3)>pitchThreshold,.5,1);
    hang.label = 'hloc';
    hang.key = 'h';
    hang.sampleRate = pch.sampleRate;    
    hang.filename = [Trial.filebase '.sst.hloc.h.mat'];
    hang.resample(locPeriods.sampleRate);
    Stc.states{end+1} = hang&wind.data;
    
end



% SPLIT paused (pause.p) state into high and low                               
sempty = isempty(Stc{'j'})||isempty(Stc{'q'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('j')) = [];
        Stc.states(Stc.gsi('q')) = [];
    end
    
    if isempty(pitchThreshold)    
% RETRIVE threshold from on pitch-rhm analysis
%         (MTA/analysis/rhythmic_head_motion/compute_session_rhm_distribution.m)
        fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
        if ~exist(fpath,'file'),
            bhv_rhm_distrb(Trial);
        end
        afig = hgload(fpath);
        rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','compute_session_rhm_distribution-hangle'));    
        rhm_distrb = rhm_distrb(1);
        figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
        mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
        mrhmp(mrhmp==0) = nan;
        rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
            +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
        pitchThreshold = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
        delete(afig);
    end
    
% INITIALIZE new states' objects with copy
    pausePeriods = Stc{'p'};
    wind = pausePeriods.copy();
    lang = pausePeriods.copy();
    hang = pausePeriods.copy();

% CREATE low pause state        
    lang.data = ThreshCross(pch(:,3)<pitchThreshold,.5,1);
    lang.label = 'lpause';
    lang.key = 'q';
    lang.sampleRate = pch.sampleRate;    
    lang.filename = [Trial.filebase '.sst.lpause.q.mat'];
    lang.resample(pausePeriods.sampleRate);
    Stc.states{end+1} = lang&wind.data;

% CREATE high pause state    
    hang.data = ThreshCross(pch(:,3)>pitchThreshold,.5,1);
    hang.label = 'hpause';
    hang.key = 'j';
    hang.sampleRate = pch.sampleRate;    
    hang.filename = [Trial.filebase '.sst.hpause.j.mat'];
    hang.resample(pausePeriods.sampleRate);    
    Stc.states{end+1} = hang&wind.data;
    
end


% END MAIN -----------------------------------------------------------------------------------------