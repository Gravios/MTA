function Trial = labelAuxBhv(Trial,Stc,varargin)
[Stc,angThresh,overwrite] = DefaultArgs(varargin,{[],-.45,false});

if isempty(Stc),
    Stc = Trial.stc.copy;
end


sempty = isempty(Stc{'t'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('t')) = [];
    end
    Stc.states{end+1} = theta(Trial);
end

sempty = isempty(Stc{'p'})||isempty(Stc{'c'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('p')) = [];
        Stc.states(Stc.gsi('c')) = [];
    end
    
    atype = 'bhv_rhm_distrb';
    ftype = 'RHM_psd_distrib_height_hangle';
    fpath = fullfile(Trial.path.data,'figures',atype,ftype,[ ftype,'-' Trial.filebase '.fig']);
    if ~exist(fpath,'file'),
        bhv_rhm_distrb(Trial);
    end
    afig = hgload(fpath);
    rhm_distrb = get(findobj(findall(get(afig,'children')),'tag','bhv_rhm_distrb-hangle'));
    rhm_distrb = get(rhm_distrb.Children);
    %figure,plot(rhm_distrb.XData,nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1)));    
    mrhmp = nanmean(rhm_distrb.CData(rhm_distrb.YData>6&rhm_distrb.YData<13,:,1));
    mrhmp(mrhmp==0) = nan;
    angThresh = rhm_distrb.XData(find(mrhmp<.5,1,'first'));
    
    % Split walking state into low and high walk
    Stc.states = cat(2,Stc.states,walk_ang(Trial,angThresh));
end


Stc.save(1);
Trial.stc = Stc;
Trial.save;


 

