function Trial = labelAuxBhv(Trial,varargin)
[Stc,angThresh,dthresh,overwrite] = DefaultArgs(varargin,{[],-.45,0.1,false});

if isempty(Stc),
    Stc = Trial.stc.copy;
end

% High and low walk
sempty = isempty(Stc{'h'})||isempty(Stc{'l'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('h')) = [];
        Stc.states(Stc.gsi('l')) = [];
    end
    
    fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
    if ~exist(fpath,'file'),
        bhv_rhm_distrb(Trial);
    end

    afig = hgload(fpath);
    rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','bhv_rhm_distrb-hangle'));
    rhm_distrb = rhm_distrb(1);
    figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
    mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
    mrhmp(mrhmp==0) = nan;
    rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
        +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
    angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
    delete(afig);
    % Split walking state into low and high walk
    wind = Stc{'w'}.copy;
    lang = Stc{'w'}.copy;
    hang = Stc{'w'}.copy;    

    ang = create(MTADang,Trial,Trial.load('xyz'));

    lang.data = ThreshCross(ang(:,5,7,2)<angThresh,.5,round(dthresh/ang.sampleRate));
    lang.label = 'lwalk';
    lang.key = 'l';
    lang.filename = [Trial.filebase '.sst.lwalk.l.mat'];
    Stc.states{end+1} = lang&wind.data;

    hang.data = ThreshCross(ang(:,5,7,2)>angThresh,.5,round(dthresh/ang.sampleRate));
    hang.label = 'hwalk';
    hang.key = 'h';
    hang.filename = [Trial.filebase '.sst.hwalk.h.mat'];
    Stc.states{end+1} = hang&wind.data;
    
end

% High and low turn
sempty = isempty(Stc{'o'})||isempty(Stc{'u'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('o')) = [];
        Stc.states(Stc.gsi('u')) = [];
    end
    
    fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
    if ~exist(fpath,'file'),
        bhv_rhm_distrb(Trial);
    end

    afig = hgload(fpath);
    rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','bhv_rhm_distrb-hangle'));
    rhm_distrb = rhm_distrb(1);
    figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
    mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
    mrhmp(mrhmp==0) = nan;
    rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
        +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
    angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
    delete(afig);
    % Split walking state into low and high walk
    wind = Stc{'t'}.copy;
    lang = Stc{'t'}.copy;
    hang = Stc{'t'}.copy;    

    ang = create(MTADang,Trial,Trial.load('xyz'));

    lang.data = ThreshCross(ang(:,5,7,2)<angThresh,.5,round(dthresh/ang.sampleRate));
    lang.label = 'lturn';
    lang.key = 'u';
    lang.filename = [Trial.filebase '.sst.lturn.u.mat'];
    Stc.states{end+1} = lang&wind.data;

    hang.data = ThreshCross(ang(:,5,7,2)>angThresh,.5,round(dthresh/ang.sampleRate));
    hang.label = 'hturn';
    hang.key = 'o';
    hang.filename = [Trial.filebase '.sst.hturn.o.mat'];
    Stc.states{end+1} = hang&wind.data;
    
end

% High and low pause 
sempty = isempty(Stc{'j'})||isempty(Stc{'q'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('j')) = [];
        Stc.states(Stc.gsi('q')) = [];
    end
    
    fpath = fullfile(Trial.spath,'figures','RHM_psd_distrib_height_hangle.fig');
    if ~exist(fpath,'file'),
        bhv_rhm_distrb(Trial);
    end

    afig = hgload(fpath);
    rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','bhv_rhm_distrb-hangle'));
    rhm_distrb = rhm_distrb(1);
    figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));   
    mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
    mrhmp(mrhmp==0) = nan;
    rhmThresh = nanmean(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)))...
        +0.5*nanstd(nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));
    angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
    delete(afig);
    % Split walking state into low and high walk
    wind = Stc{'p'}.copy;
    lang = Stc{'p'}.copy;
    hang = Stc{'p'}.copy;    

    ang = create(MTADang,Trial,Trial.load('xyz'));

    lang.data = ThreshCross(ang(:,5,7,2)<angThresh,.5,round(dthresh/ang.sampleRate));
    lang.label = 'lpause';
    lang.key = 'q';
    lang.filename = [Trial.filebase '.sst.lpause.q.mat'];
    Stc.states{end+1} = lang&wind.data;

    hang.data = ThreshCross(ang(:,5,7,2)>angThresh,.5,round(dthresh/ang.sampleRate));
    hang.label = 'hpause';
    hang.key = 'j';
    hang.filename = [Trial.filebase '.sst.hpause.j.mat'];
    Stc.states{end+1} = hang&wind.data;
    
end


Stc.save(1);
Trial.stc = Stc;
Trial.save;

