function Trial = label_aux_bhv_reduced(Trial,varargin)
[Stc,angThresh,dthresh,overwrite] = DefaultArgs(varargin,{[],[],0.1,false});

if isempty(Stc),
    Stc = Trial.stc.copy;
elseif ischar(Stc),
    Stc = Trial.load('stc',Stc);
end

% High and low walk
sempty = isempty(Stc{'h'})||isempty(Stc{'l'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('h')) = [];
        Stc.states(Stc.gsi('l')) = [];
    end
    

    if isempty(angThresh)
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
        angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
        delete(afig);
    end
    
    
    % Split walking state into low and high walk
    wind = Stc{'x'}.copy;
    lang = Stc{'x'}.copy;
    hang = Stc{'x'}.copy;    

    xyz = Trial.load('xyz');
    xyz.filter('ButFilter',3,1.4,'low');    
    ang = create(MTADang,Trial,xyz);

    lang.data = ThreshCross(ang(:,5,7,2)<angThresh,.5,round(dthresh/ang.sampleRate));
    lang.label = 'lloc';
    lang.key = 'l';
    lang.filename = [Trial.filebase '.sst.lloc.l.mat'];
    Stc.states{end+1} = lang&wind.data;

    hang.data = ThreshCross(ang(:,5,7,2)>angThresh,.5,round(dthresh/ang.sampleRate));
    hang.label = 'hloc';
    hang.key = 'h';
    hang.filename = [Trial.filebase '.sst.hloc.h.mat'];
    Stc.states{end+1} = hang&wind.data;
    
end

% High and low pause 
sempty = isempty(Stc{'j'})||isempty(Stc{'q'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('j')) = [];
        Stc.states(Stc.gsi('q')) = [];
    end
    
    if isempty(angThresh)    
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
        angThresh = rhm_distrb.XData(find(mrhmp<rhmThresh,1,'first'));
        delete(afig);
    end
    
    % Split walking state into low and high walk
    wind = Stc{'p'}.copy;
    lang = Stc{'p'}.copy;
    hang = Stc{'p'}.copy;    

    xyz = Trial.load('xyz');
    xyz.filter('ButFilter',3,1.4,'low');        
    ang = create(MTADang,Trial,xyz);

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

